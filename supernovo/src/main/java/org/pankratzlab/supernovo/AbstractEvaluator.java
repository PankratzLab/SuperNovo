package org.pankratzlab.supernovo;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import org.pankratzlab.supernovo.output.DeNovoResult;
import org.pankratzlab.supernovo.pileup.Depth;
import org.pankratzlab.supernovo.pileup.Pileup;
import org.pankratzlab.supernovo.pileup.PileupCache;
import com.google.common.base.Optional;
import com.google.common.base.Predicates;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.ConcurrentHashMultiset;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.LinkedHashMultiset;
import com.google.common.collect.MoreCollectors;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public abstract class AbstractEvaluator implements Evaluator {

  protected static final String SER_EXTENSION = ".SuperNovoResultList.ser.gz";
  protected static final String VCF_EXTENSION = ".SuperNovoResults.vcf.gz";

  private final Multiset<String> contigLogCount = ConcurrentHashMultiset.create();
  private final Function<File, VCFHeader> vcfHeaderCache =
      CacheBuilder.newBuilder()
          .softValues()
          .build(
              CacheLoader.from(
                  vcf -> {
                    try (VCFFileReader vcfHeaderReader = new VCFFileReader(vcf)) {
                      return vcfHeaderReader.getFileHeader();
                    }
                  }));

  private ConcurrentHashMap<ReferencePosition, Optional<DeNovoResult>> results;
  private final PileupCache childPileups;

  private void maybeLogParseProgress(VariantContext vc) {
    int processedVariants = contigLogCount.count(vc.getContig());
    if (processedVariants % 10000 == 0 && processedVariants != 0) {
      App.LOG.info(
          "Parsed "
              + contigLogCount.count(vc.getContig())
              + " positions on contig "
              + vc.getContig());
    }
    contigLogCount.add(vc.getContig());
  }

  private void serializeChunkedProgress(File chunkedSerOutput, Thread generateResultsThread) {
    long lastSerialize = System.currentTimeMillis();
    while (generateResultsThread.isAlive()) {
      if (System.currentTimeMillis() - lastSerialize > 600000) {
        File tempChunkedSerOutput = new File(chunkedSerOutput.getPath() + "_TEMP");
        serializeResults(tempChunkedSerOutput);
        if (!tempChunkedSerOutput.renameTo(chunkedSerOutput)) {
          App.LOG.error("Failed to overwrite temp chunked output, chunking may not be reloadable");
        }
        lastSerialize = System.currentTimeMillis();
      } else {
        try {
          Thread.sleep(60000);
        } catch (InterruptedException e) {
          Thread.currentThread().interrupt();
        }
      }
    }
  }

  private void generateResults(File vcf) {
    LoadingCache<Thread, VCFFileReader> perThreadReaders =
        CacheBuilder.newBuilder().build(CacheLoader.from(t -> new VCFFileReader(vcf)));
    Supplier<VCFFileReader> getVCFReader =
        () -> perThreadReaders.getUnchecked(Thread.currentThread());
    long time = System.currentTimeMillis();
    App.LOG.info("Parsing variants from vcf");
    ImmutableSet<ReferencePosition> variantsToInclude =
        genomeBins(vcf)
            .parallel()
            .map(loc -> getVCFReader.get().query(loc).stream())
            .flatMap(s -> s.peek(this::maybeLogParseProgress).filter(this::keepVariant))
            .map(this::generatePosition)
            .filter(Optional::isPresent)
            .map(Optional::get)
            .collect(ImmutableSet.toImmutableSet());
    if (!results.isEmpty()) {
      int prevResults = results.size();
      results.keySet().retainAll(variantsToInclude);
      App.LOG.info(
          "Dropped "
              + (prevResults - results.size())
              + " results (of "
              + prevResults
              + " total) from previously computed that are no longer retained from vcf");
    }
    ImmutableSet<ReferencePosition> variantsRemaining =
        variantsToInclude
            .stream()
            .filter(Predicates.not(results::containsKey))
            .collect(ImmutableSet.toImmutableSet());
    App.LOG.info("Parsed variants in " + (System.currentTimeMillis() - time) / 1000 + " seconds");
    perThreadReaders.asMap().values().forEach(VCFFileReader::close);

    App.LOG.info("Evaluating " + variantsRemaining.size() + " variants for de novo mutations");
    if (!results.isEmpty()) {
      App.LOG.info(results.size() + " variants previously evaluated");
    }
    Thread evaluateVariantsThreads =
        new Thread(
            () ->
                variantsRemaining
                    .stream()
                    .parallel()
                    .map(r -> ImmutableMap.of(r, evaluate(r)))
                    .forEach(results::putAll));
    evaluateVariantsThreads.start();
    new Thread(() -> logProgress(variantsToInclude.size(), evaluateVariantsThreads)).start();
    try {
      evaluateVariantsThreads.join();
    } catch (InterruptedException e) {
      Thread.currentThread().interrupt();
    }
    App.LOG.info(
        "Finished evaluating "
            + results.size()
            + " (of "
            + variantsToInclude.size()
            + " expected) variants for de novo mutations");
  }

  private void logProgress(int totalVariantsCount, Thread evaluateVariantsThreads) {
    long lastProgressLog = 0L;
    while (evaluateVariantsThreads.isAlive()) {
      if (System.currentTimeMillis() - lastProgressLog > 600000) {
        App.LOG.info(
            "Processed " + results.size() + " variants (of " + totalVariantsCount + " total");
        lastProgressLog = System.currentTimeMillis();
      }
      try {
        Thread.sleep(60000);
      } catch (InterruptedException e) {
        Thread.currentThread().interrupt();
      }
    }
  }

  public AbstractEvaluator(File bam) {
    super();
    this.childPileups = new PileupCache(bam);
  }

  @Override
  public void run(File vcf, File output) throws IOException, ClassNotFoundException {
    File vcfOutput = formVCFOutput(output);
    File serOutput = formSerializedOutput(output);
    File chunkedSerOutput = formChunkedSerializedOutput(output);
    Optional<ConcurrentHashMap<ReferencePosition, Optional<DeNovoResult>>> prevResults =
        Optional.absent();
    if (serOutput.exists()) {
      App.LOG.info("Previous serialized output already exists, loading...");
      try {
        prevResults = Optional.of(deserializeResults(serOutput));
        App.LOG.info("Serialized output loaded");
      } catch (Exception e) {
        App.LOG.error("Error loading serialized results, regenerating", e);
      }
    }
    if (chunkedSerOutput.exists() && !prevResults.isPresent()) {
      App.LOG.info("Previous chunked progress serialized output already exists, loading...");
      try {
        prevResults = Optional.of(deserializeResults(chunkedSerOutput));
        App.LOG.info("Serialized output loaded");
      } catch (Exception e) {
        App.LOG.error("Error loading serialized results, regenerating", e);
      }
    }
    results = prevResults.or(ConcurrentHashMap::new);
    Thread generateResultsThread = new Thread(() -> this.generateResults(vcf));
    generateResultsThread.start();
    new Thread(() -> this.serializeChunkedProgress(chunkedSerOutput, generateResultsThread))
        .start();
    try {
      generateResultsThread.join();
    } catch (InterruptedException e) {
      Thread.currentThread().interrupt();
    }
    ImmutableList<DeNovoResult> resultsList =
        results
            .values()
            .stream()
            .filter(Optional::isPresent)
            .map(Optional::get)
            .sorted(Comparator.comparing(DeNovoResult::getPos))
            .collect(ImmutableList.toImmutableList());
    DeNovoResult.retrieveAnnos(
        resultsList, vcfOutput, vcfHeaderCache.apply(vcf).getSequenceDictionary());
    serializeResults(serOutput);
    summarizeResults(resultsList, formSummarizedOutput(output));
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(output)))) {
      resultsList.stream().findFirst().map(DeNovoResult::generateHeader).ifPresent(writer::println);
      resultsList.stream().map(DeNovoResult::generateLine).forEachOrdered(writer::println);
    }
  }

  protected abstract DeNovoResult generateDeNovoResult(
      ReferencePosition pos, Pileup childPile, PileupCache childPileups);

  protected boolean keepVariant(VariantContext vc) {
    Genotype geno = vc.getGenotype(App.getInstance().getChildID());
    return isSingleNonRef(geno);
  }

  protected Allele getPossibleDnAlt(VariantContext vc) {
    Genotype geno = vc.getGenotype(App.getInstance().getChildID());
    return geno.getAlleles()
        .stream()
        .filter(Predicates.not(vc.getReference()::equals))
        .collect(MoreCollectors.onlyElement());
  }

  protected DeNovoResult.Sample generateSample(
      String id, ReferencePosition pos, Pileup pileup, Pileup childPile) {
    return new DeNovoResult.Sample(
        id, pileup, pos, childPile.getDepth().getA1(), childPile.getDepth().getA2());
  }

  private Stream<Locatable> genomeBins(File vcf) {
    final int binSize = 100000;
    return vcfHeaderCache
        .apply(vcf)
        .getContigLines()
        .stream()
        .map(VCFContigHeaderLine::getSAMSequenceRecord)
        .flatMap(
            r ->
                IntStream.rangeClosed(0, r.getSequenceLength() / binSize)
                    .mapToObj(
                        i ->
                            new Interval(r.getSequenceName(), i * binSize + 1, (i + 1) * binSize)));
  }

  private void serializeResults(File output) {
    results.size();
    try (ObjectOutputStream oos =
        new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(output)))) {
      oos.writeObject(results);
    } catch (IOException e) {
      App.LOG.error(e);
    }
  }

  private void summarizeResults(ImmutableList<DeNovoResult> results, File output)
      throws IOException {
    Multiset<String> counts = LinkedHashMultiset.create();
    for (DeNovoResult dnr : results) {
      if (dnr.superNovo) {
        counts.add("supernovo");
        counts.add(dnr.snpeffGene + "_AnyImpact");
        counts.add(dnr.snpeffImpact);
        if ("MODERATE".equals(dnr.snpeffImpact) || "HIGH".equals(dnr.snpeffImpact)) {
          counts.add("supernovo_damaging");
          counts.add(dnr.snpeffGene);
          if (!dnr.dnIsRef.or(Boolean.FALSE)) counts.add("supernovo_damaging_nonref");
        }
      }
    }
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(output)))) {
      counts
          .entrySet()
          .stream()
          .map(e -> e.getElement() + "\t" + e.getCount())
          .forEach(writer::println);
    }
  }

  private boolean isSingleNonRef(Genotype geno) {
    return (geno.getPloidy() == 1 || (geno.isHet() && !geno.isHetNonRef()))
        && geno.getAlleles().stream().mapToInt(Allele::length).allMatch(i -> i == 1);
  }

  private Optional<ReferencePosition> generatePosition(VariantContext vc) {
    Allele ref = vc.getReference();
    Allele alt = getPossibleDnAlt(vc);
    try {
      return Optional.of(ReferencePosition.fromVariantContext(vc, ref, alt));
    } catch (IllegalArgumentException iae) {
      App.LOG.error("Failed to generate ReferencePosition for variant " + vc, iae);
      return Optional.absent();
    }
  }

  private Optional<DeNovoResult> evaluate(ReferencePosition pos) {
    Pileup childPile = childPileups.get(pos);
    if (looksVariant(childPile.getDepth())) {
      return Optional.of(generateDeNovoResult(pos, childPile, childPileups));
    }
    return Optional.absent();
  }

  private static File formSerializedOutput(File textOutput) {
    String path = textOutput.getPath();
    return new File(path.substring(0, path.lastIndexOf('.')) + SER_EXTENSION);
  }

  private static File formChunkedSerializedOutput(File textOutput) {
    return new File(formSerializedOutput(textOutput).getPath() + "_CHUNKED");
  }

  private static File formVCFOutput(File textOutput) {
    String path = textOutput.getPath();
    return new File(path.substring(0, path.lastIndexOf('.')) + VCF_EXTENSION);
  }

  private static File formSummarizedOutput(File textOutput) {
    String path = textOutput.getPath();
    return new File(path.substring(0, path.lastIndexOf('.')) + ".summary.txt");
  }

  private static Set<PileAllele> possibleAlleles(Pileup pileup) {
    return pileup
        .getBaseCounts()
        .entrySet()
        .stream()
        .filter(
            e ->
                e.getCount() > App.getInstance().getMaxMiscallWeight()
                    || pileup.getBaseFractions().get(e.getElement())
                        > App.getInstance().getMaxMiscallFrac())
        .map(Multiset.Entry::getElement)
        .collect(ImmutableSet.toImmutableSet());
  }

  public static boolean looksBiallelic(Pileup pileup) {
    return looksVariant(pileup.getDepth()) && !moreThanTwoViableAlleles(pileup);
  }

  public static boolean looksVariant(Depth depth) {
    return depth.getBiAlleles().size() == 2
        && depth.weightedBiallelicDepth() >= App.getInstance().getMinDepth()
        && passesAllelicFrac(depth)
        && passesAllelicDepth(depth, App.getInstance().getMinAllelicDepth());
  }

  public static boolean passesAllelicFrac(Depth depth) {
    return depth.weightedMinorAlleleFraction() >= App.getInstance().getMinAllelicFrac();
  }

  public static boolean passesAllelicDepth(Depth depth, double minDepth) {
    return Arrays.stream(Depth.Allele.values())
        .mapToDouble(depth::allelicRawDepth)
        .allMatch(d -> d >= minDepth);
  }

  public static boolean moreThanTwoViableAlleles(Pileup pileup) {
    return possibleAlleles(pileup).size() > 2;
  }

  public static boolean looksDenovo(
      Pileup childPileup, Optional<Pileup> p1Pileup, Optional<Pileup> p2Pileup) {
    return dnAllele(childPileup, p1Pileup, p2Pileup).isPresent();
  }

  public static Optional<PileAllele> dnAllele(
      Pileup childPileup, Optional<Pileup> p1Pileup, Optional<Pileup> p2Pileup) {
    Set<PileAllele> parentalAlleles =
        Stream.of(p1Pileup, p2Pileup)
            .filter(Optional::isPresent)
            .map(Optional::get)
            .map(pileup -> possibleAlleles(pileup))
            .flatMap(Set::stream)
            .collect(ImmutableSet.toImmutableSet());
    try {
      return Optional.fromJavaUtil(
          Sets.difference(childPileup.getDepth().getBiAlleles(), parentalAlleles)
              .stream()
              .collect(MoreCollectors.toOptional()));
    } catch (IllegalArgumentException iae) {
      if (p1Pileup.isPresent() && p2Pileup.isPresent()) {
        App.LOG.warn(
            "Multiple alleles at site appear De Novo for child: "
                + childPileup
                + ", with parents: "
                + p1Pileup
                + " and "
                + p2Pileup);
      }
      return Optional.absent();
    }
  }

  @SuppressWarnings("unchecked")
  public static ConcurrentHashMap<ReferencePosition, Optional<DeNovoResult>> deserializeResults(
      File input) throws IOException, ClassNotFoundException {
    try (ObjectInputStream ois =
        new ObjectInputStream(
            new BufferedInputStream(new GZIPInputStream(new FileInputStream(input))))) {
      return (ConcurrentHashMap<ReferencePosition, Optional<DeNovoResult>>) ois.readObject();
    }
  }
}
