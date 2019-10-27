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
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import org.apache.logging.log4j.LogManager;
import org.pankratzlab.supernovo.output.DeNovoResult;
import org.pankratzlab.supernovo.output.OutputFields;
import org.pankratzlab.supernovo.pileup.Depth;
import org.pankratzlab.supernovo.pileup.Pileup;
import org.pankratzlab.supernovo.pileup.SAMPositionQueryOverlap;
import com.google.common.base.Optional;
import com.google.common.base.Predicates;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.MoreCollectors;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFileReader;

public class TrioEvaluator {

  private static final int READ_LENGTH = 150;
  private static final int MIN_DEPTH = 10;
  private static final int MIN_ALLELIC_DEPTH = 4;
  private static final double MIN_ALLELIC_FRAC = 0.1;
  private static final double MAX_MISCALL_RATIO = 0.05;
  private static final double MAX_MISCALL_WEIGHT = 1.0;
  private static final CacheBuilder<Object, Object> PILEUP_CACHE_BUILDER =
      CacheBuilder.newBuilder().maximumSize(READ_LENGTH * 2L);

  private final String childID;
  private final String parent1ID;
  private final String parent2ID;

  private final LoadingCache<GenomePosition, Pileup> childPileups;
  private final LoadingCache<GenomePosition, Pileup> p1Pileups;
  private final LoadingCache<GenomePosition, Pileup> p2Pileups;

  /**
   * @param childBam {@link SamReader} of child to evluate for de novo variants
   * @param parent1Bam {@link SamReader} of one parent for child
   * @param parent2Bam {@link SamReader} of second parent for child
   */
  public TrioEvaluator(
      File childBam,
      String childID,
      File parent1Bam,
      String parent1ID,
      File parent2Bam,
      String parent2ID) {
    super();
    this.childID = childID;
    this.parent1ID = parent1ID;
    this.parent2ID = parent2ID;

    this.childPileups = PILEUP_CACHE_BUILDER.build(queryingPileupLoader(childBam));
    this.p1Pileups = PILEUP_CACHE_BUILDER.build(queryingPileupLoader(parent1Bam));
    this.p2Pileups = PILEUP_CACHE_BUILDER.build(queryingPileupLoader(parent2Bam));
  }

  private static CacheLoader<GenomePosition, Pileup> queryingPileupLoader(final File bam) {
    return new CacheLoader<GenomePosition, Pileup>() {

      @Override
      public Pileup load(GenomePosition pos) {
        try (SAMPositionQueryOverlap spqo = new SAMPositionQueryOverlap(bam, pos)) {
          return new Pileup(spqo.getRecords(), pos);
        }
      }
    };
  }

  public void reportDeNovos(File vcf, File output) throws IOException {
    File serOutput = formSerializedOutput(output);
    final ImmutableList<DeNovoResult> results;
    if (serOutput.exists()) {
      App.LOG.info("Serialized output already exists, loading...");
      try {
        results = deserializeResults(serOutput);
        App.LOG.info("Serialized output loaded");
      } catch (ClassNotFoundException e) {
        App.LOG.error("Error loading serialized results", e);
        return;
      }
    } else {
      results =
          perContigVCFReaders(vcf)
              .parallel()
              .flatMap(
                  s ->
                      s.filter(this::keepVariant)
                          .map(this::generatePosition)
                          .filter(Optional::isPresent)
                          .map(Optional::get)
                          .map(this::evaluate)
                          .filter(Optional::isPresent)
                          .map(Optional::get))
              .collect(ImmutableList.toImmutableList());
    }
    DeNovoResult.retrieveAnnos(results);
    serializeResults(results, serOutput);
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(output)))) {
      writer.println(OutputFields.generateHeader(DeNovoResult.class));
      results.stream().map(DeNovoResult::generateLine).forEachOrdered(writer::println);
    }
  }

  private Stream<Stream<VariantContext>> perContigVCFReaders(File vcf) {
    return new VCFFileReader(vcf)
        .getFileHeader()
        .getContigLines()
        .stream()
        .map(VCFContigHeaderLine::getID)
        .map(v -> new VCFFileReader(vcf).query(v, 1, Integer.MAX_VALUE))
        .map(CloseableIterator::stream);
  }

  private static File formSerializedOutput(File textOutput) {
    String path = textOutput.getPath();
    return new File(path.substring(0, path.lastIndexOf('.')) + ".DeNovoResultList.ser.gz");
  }

  private void serializeResults(ImmutableList<DeNovoResult> results, File output) {
    try (ObjectOutputStream oos =
        new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(output)))) {
      oos.writeObject(results);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  @SuppressWarnings("unchecked")
  private ImmutableList<DeNovoResult> deserializeResults(File input)
      throws IOException, ClassNotFoundException {
    try (ObjectInputStream ois =
        new ObjectInputStream(
            new BufferedInputStream(new GZIPInputStream(new FileInputStream(input))))) {
      return (ImmutableList<DeNovoResult>) ois.readObject();
    }
  }

  private boolean keepVariant(VariantContext vc) {
    Genotype geno = vc.getGenotype(childID);
    return geno.isHet()
        && !geno.isHetNonRef()
        && geno.getAlleles().stream().mapToInt(Allele::length).anyMatch(i -> i == 1);
  }

  private Optional<ReferencePosition> generatePosition(VariantContext vc) {
    Allele ref = vc.getReference();
    Genotype geno = vc.getGenotype(childID);
    Allele alt =
        geno.getAlleles()
            .stream()
            .filter(Predicates.not(vc.getReference()::equals))
            .collect(MoreCollectors.onlyElement());
    try {
      return Optional.of(ReferencePosition.fromVariantContext(vc, ref, alt));
    } catch (IllegalArgumentException iae) {
      LogManager.getLogger(App.class)
          .error("Failed to generate ReferencePosition for variant", iae);
      return Optional.absent();
    }
  }

  private Optional<DeNovoResult> evaluate(ReferencePosition pos) {
    Pileup childPile = childPileups.getUnchecked(pos);
    if (looksVariant(childPile.getDepth())) {
      return Optional.of(
          new DeNovoResult(
              pos,
              new HaplotypeEvaluator(
                      childPile,
                      childPileups::getUnchecked,
                      p1Pileups::getUnchecked,
                      p2Pileups::getUnchecked)
                  .haplotypeConcordance(),
              generateSample(childID, pos, childPile, childPile),
              generateSample(parent1ID, pos, p1Pileups.getUnchecked(pos), childPile),
              generateSample(parent2ID, pos, p2Pileups.getUnchecked(pos), childPile)));
    }
    return Optional.absent();
  }

  public static boolean looksBiallelic(Pileup pileup) {
    return looksVariant(pileup.getDepth()) && !moreThanTwoViableAlleles(pileup);
  }

  public static boolean looksVariant(Depth depth) {
    return depth.getBiAlleles().size() == 2
        && depth.weightedBiallelicDepth() >= MIN_DEPTH
        && depth.weightedMinorAlleleFraction() >= MIN_ALLELIC_FRAC
        && Arrays.stream(Depth.Allele.values())
            .mapToDouble(depth::allelicWeightedDepth)
            .allMatch(d -> d >= MIN_ALLELIC_DEPTH);
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
    List<Pileup> parentPileups = ImmutableList.of(p1Pileup, p2Pileup);
    Set<PileAllele> parentalAlleles =
        parentPileups
            .stream()
            .map(TrioEvaluator::possibleAlleles)
            .flatMap(Set::stream)
            .collect(ImmutableSet.toImmutableSet());
    return !Sets.difference(childPileup.getDepth().getBiAlleles(), parentalAlleles).isEmpty();
  }
}
