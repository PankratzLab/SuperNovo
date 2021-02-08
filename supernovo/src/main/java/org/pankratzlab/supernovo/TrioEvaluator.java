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
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import org.apache.logging.log4j.LogManager;
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
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class TrioEvaluator {

  protected static final String SER_EXTENSION = ".DeNovoResultList.ser.gz";
  protected static final String VCF_EXTENSION = ".DeNovoResults.vcf.gz";

  private static final int VCF_DN_MAX_WEIGHT = 4;
  private static final int MIN_DEPTH = 10;
  private static final int MIN_ALLELIC_DEPTH = 4;
  private static final double MIN_ALLELIC_FRAC = 0.1;
  private static final double MAX_MISCALL_RATIO = 0.05;
  private static final double MAX_MISCALL_WEIGHT = 1.0;

  private final String childID;
  private final String parent1ID;
  private final String parent2ID;
  private final String snpEffGenome;

  private final PileupCache childPileups;
  private final PileupCache p1Pileups;
  private final PileupCache p2Pileups;

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

  private ConcurrentHashMap<ReferencePosition, Optional<DeNovoResult>> deNovoResults;

  /**
   * @param childBam {@link SamReader} of child to evluate for de novo variants
   * @param parent1Bam {@link SamReader} of one parent for child
   * @param parent2Bam {@link SamReader} of second parent for child
   * @param snpEffGenome genome build argument to supply SnpEff
   */
  public TrioEvaluator(
      File childBam,
      String childID,
      File parent1Bam,
      String parent1ID,
      File parent2Bam,
      String parent2ID,
      String snpEffGenome) {
    super();
    this.childID = childID;
    this.parent1ID = parent1ID;
    this.parent2ID = parent2ID;
    this.snpEffGenome = snpEffGenome;

    this.childPileups = new PileupCache(childBam);
    this.p1Pileups = new PileupCache(parent1Bam);
    this.p2Pileups = new PileupCache(parent2Bam);
  }

  private void generateResults(File vcf) {
    LoadingCache<Thread, VCFFileReader> perThreadReaders =
        CacheBuilder.newBuilder().build(CacheLoader.from(t -> new VCFFileReader(vcf)));
    Supplier<VCFFileReader> getVCFReader =
        () -> perThreadReaders.getUnchecked(Thread.currentThread());
    long time = System.currentTimeMillis();
    App.LOG.info("Parsing variants from gvcf");
    ImmutableSet<ReferencePosition> variantsToInclude =
        genomeBins(vcf)
            .parallel()
            .map(loc -> getVCFReader.get().query(loc).stream())
            .flatMap(s -> s.map(this::maybeLogParseProgress).filter(this::keepVariant))
            .map(this::generatePosition)
            .filter(Optional::isPresent)
            .map(Optional::get)
            .collect(ImmutableSet.toImmutableSet());
    int prevResults = deNovoResults.size();
    deNovoResults.keySet().retainAll(variantsToInclude);
    App.LOG.info(
        "Dropped "
            + (prevResults - deNovoResults.size())
            + " results (of "
            + prevResults
            + " total) from previously computed that are no longer retained from vcf");
    ImmutableSet<ReferencePosition> variantsRemaining =
        variantsToInclude
            .stream()
            .filter(Predicates.not(deNovoResults::containsKey))
            .collect(ImmutableSet.toImmutableSet());
    App.LOG.info("Parsed variants in " + (System.currentTimeMillis() - time) + " seconds");
    perThreadReaders.asMap().values().forEach(VCFFileReader::close);

    App.LOG.info("Evaluating " + variantsRemaining.size() + " variants for de novo mutations");
    if (!deNovoResults.isEmpty()) {
      App.LOG.info(deNovoResults.size() + " variants previously evaluated");
    }
    Thread evaluateVariantsThreads =
        new Thread(
            () ->
                variantsRemaining
                    .stream()
                    .parallel()
                    .map(r -> ImmutableMap.of(r, evaluate(r)))
                    .forEach(deNovoResults::putAll));
    evaluateVariantsThreads.start();
    long lastProgressLog = 0L;
    while (evaluateVariantsThreads.isAlive()) {
      if (System.currentTimeMillis() - lastProgressLog > 600000) {
        App.LOG.info(
            "Processed "
                + deNovoResults.size()
                + " variants (of "
                + variantsRemaining.size()
                + " total");
      }
      try {
        Thread.sleep(60000);
      } catch (InterruptedException e) {
        Thread.currentThread().interrupt();
      }
    }

    App.LOG.info(
        "Finished evaluating " + variantsRemaining.size() + " variants for de novo mutations");
  }

  private VariantContext maybeLogParseProgress(VariantContext vc) {
    if (contigLogCount.count(vc.getContig()) % 10000 == 0) {
      App.LOG.info(
          "Parsed "
              + contigLogCount.count(vc.getContig())
              + " positions on contig "
              + vc.getContig());
    }
    contigLogCount.add(vc.getContig());
    return vc;
  }

  public void reportDeNovos(File vcf, File output) throws IOException, ClassNotFoundException {
    File vcfOutput = formVCFOutput(output);
    File serOutput = formSerializedOutput(output);
    File chunkedSerOutput = formChunkedSerializedOutput(output);
    deNovoResults = new ConcurrentHashMap<>();
    if (serOutput.exists()) {
      App.LOG.info("Previous serialized output already exists, loading...");
      try {
        deNovoResults = deserializeResults(serOutput);
        App.LOG.info("Serialized output loaded");
      } catch (Exception e) {
        App.LOG.error("Error loading serialized results, regenerating", e);
      }
    }
    if (chunkedSerOutput.exists()) {
      App.LOG.info("Previous chunked progress serialized output already exists, loading...");
      try {
        deNovoResults = deserializeResults(chunkedSerOutput);
        App.LOG.info("Serialized output loaded");
      } catch (Exception e) {
        App.LOG.error("Error loading serialized results, regenerating", e);
      }
    }
    Thread generateResultsThread = new Thread(() -> this.generateResults(vcf));
    generateResultsThread.start();
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
    ImmutableList<DeNovoResult> resultsList =
        deNovoResults
            .values()
            .stream()
            .filter(Optional::isPresent)
            .map(Optional::get)
            .sorted(Comparator.comparing(DeNovoResult::getPos))
            .collect(ImmutableList.toImmutableList());
    DeNovoResult.retrieveAnnos(
        resultsList, vcfOutput, vcfHeaderCache.apply(vcf).getSequenceDictionary(), snpEffGenome);
    serializeResults(serOutput);
    summarizeResults(resultsList, formSummarizedOutput(output));
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(output)))) {
      resultsList.stream().findFirst().map(DeNovoResult::generateHeader).ifPresent(writer::println);
      resultsList.stream().map(DeNovoResult::generateLine).forEachOrdered(writer::println);
    }
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

  private void serializeResults(File output) {
    deNovoResults.size();
    try (ObjectOutputStream oos =
        new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(output)))) {
      oos.writeObject(deNovoResults);
    } catch (IOException e) {
      e.printStackTrace();
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

  @SuppressWarnings("unchecked")
  protected static ConcurrentHashMap<ReferencePosition, Optional<DeNovoResult>> deserializeResults(
      File input) throws IOException, ClassNotFoundException {
    try (ObjectInputStream ois =
        new ObjectInputStream(
            new BufferedInputStream(new GZIPInputStream(new FileInputStream(input))))) {
      return (ConcurrentHashMap<ReferencePosition, Optional<DeNovoResult>>) ois.readObject();
    }
  }

  private boolean keepVariant(VariantContext vc) {
    Genotype geno = vc.getGenotype(childID);
    return isSingleNonRef(geno) && !seenInParentVCF(vc);
  }

  private boolean isSingleNonRef(Genotype geno) {
    return (geno.getPloidy() == 1 || (geno.isHet() && !geno.isHetNonRef()))
        && geno.getAlleles().stream().mapToInt(Allele::length).allMatch(i -> i == 1);
  }

  private boolean seenInParentVCF(VariantContext vc) {
    Allele altAllele = getPossibleDnAlt(vc);
    for (String parentID : ImmutableList.of(parent1ID, parent2ID)) {
      Genotype geno = vc.getGenotype(parentID);
      int parentAlleleIndex = geno.getAlleles().indexOf(altAllele);
      if (parentAlleleIndex >= 0 && geno.getAD()[parentAlleleIndex] > VCF_DN_MAX_WEIGHT) {
        return true;
      }
    }
    return false;
  }

  private Optional<ReferencePosition> generatePosition(VariantContext vc) {
    Allele ref = vc.getReference();
    Allele alt = getPossibleDnAlt(vc);
    try {
      return Optional.of(ReferencePosition.fromVariantContext(vc, ref, alt));
    } catch (IllegalArgumentException iae) {
      LogManager.getLogger(App.class)
          .error("Failed to generate ReferencePosition for variant " + vc, iae);
      return Optional.absent();
    }
  }

  private Allele getPossibleDnAlt(VariantContext vc) {
    Genotype geno = vc.getGenotype(childID);
    return geno.getAlleles()
        .stream()
        .filter(Predicates.not(vc.getReference()::equals))
        .collect(MoreCollectors.onlyElement());
  }

  private Optional<DeNovoResult> evaluate(ReferencePosition pos) {
    Pileup childPile = childPileups.get(pos);
    if (looksVariant(childPile.getDepth())) {
      return Optional.of(
          new DeNovoResult(
              pos,
              new HaplotypeEvaluator(childPile, childPileups, p1Pileups, p2Pileups)
                  .haplotypeConcordance(),
              generateSample(childID, pos, childPile, childPile),
              generateSample(parent1ID, pos, p1Pileups.get(pos), childPile),
              generateSample(parent2ID, pos, p2Pileups.get(pos), childPile)));
    }
    return Optional.absent();
  }

  public static boolean looksBiallelic(Pileup pileup) {
    return looksVariant(pileup.getDepth()) && !moreThanTwoViableAlleles(pileup);
  }

  public static boolean looksVariant(Depth depth) {
    return depth.getBiAlleles().size() == 2
        && depth.weightedBiallelicDepth() >= MIN_DEPTH
        && passesAllelicFrac(depth)
        && passesAllelicDepth(depth, MIN_ALLELIC_DEPTH);
  }

  public static boolean passesAllelicFrac(Depth depth) {
    return depth.weightedMinorAlleleFraction() >= MIN_ALLELIC_FRAC;
  }

  public static boolean passesAllelicDepth(Depth depth, double minDepth) {
    return Arrays.stream(Depth.Allele.values())
        .mapToDouble(depth::allelicRawDepth)
        .allMatch(d -> d >= minDepth);
  }

  public static boolean moreThanTwoViableAlleles(Pileup pileup) {
    return possibleAlleles(pileup).size() > 2;
  }

  private static Set<PileAllele> possibleAlleles(Pileup pileup) {
    return pileup
        .getBaseCounts()
        .entrySet()
        .stream()
        .filter(
            e ->
                e.getCount() > MAX_MISCALL_WEIGHT
                    || pileup.getBaseFractions().get(e.getElement()) > MAX_MISCALL_RATIO)
        .map(Multiset.Entry::getElement)
        .collect(ImmutableSet.toImmutableSet());
  }

  private static DeNovoResult.Sample generateSample(
      String id, ReferencePosition pos, Pileup pileup, Pileup childPile) {
    return new DeNovoResult.Sample(
        id, pileup, pos, childPile.getDepth().getA1(), childPile.getDepth().getA2());
  }

  public static boolean looksDenovo(Pileup childPileup, Pileup p1Pileup, Pileup p2Pileup) {
    return dnAllele(childPileup, p1Pileup, p2Pileup).isPresent();
  }

  public static Optional<PileAllele> dnAllele(
      Pileup childPileup, Pileup p1Pileup, Pileup p2Pileup) {
    List<Pileup> parentPileups = ImmutableList.of(p1Pileup, p2Pileup);
    Set<PileAllele> parentalAlleles =
        parentPileups
            .stream()
            .map(TrioEvaluator::possibleAlleles)
            .flatMap(Set::stream)
            .collect(ImmutableSet.toImmutableSet());
    try {
      return Optional.fromJavaUtil(
          Sets.difference(childPileup.getDepth().getBiAlleles(), parentalAlleles)
              .stream()
              .collect(MoreCollectors.toOptional()));
    } catch (IllegalArgumentException iae) {
      App.LOG.warn(
          "Multiple alleles at site appear De Novo for child: "
              + childPileup
              + ", with parents: "
              + p1Pileup
              + " and "
              + p2Pileup);
      return Optional.absent();
    }
  }
}
