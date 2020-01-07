package org.pankratzlab.supernovo.pileup;

import java.io.Serializable;
import java.util.Comparator;
import java.util.Map;
import java.util.stream.Stream;
import org.pankratzlab.supernovo.GenomePosition;
import org.pankratzlab.supernovo.PileAllele;
import org.pankratzlab.supernovo.ReferencePosition;
import org.pankratzlab.supernovo.SNPAllele;
import com.google.common.base.Optional;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.Maps;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class Pileup implements Serializable {

  /** */
  private static final long serialVersionUID = 2L;

  public static class Builder {
    private final GenomePosition position;
    private final ImmutableList<PileAllele> queriedAlleles;
    private final ImmutableSetMultimap.Builder<PileAllele, Integer> basePilesBuilder;
    private final Map<PileAllele, Double> weightedDepth;
    private final ImmutableMultiset.Builder<PileAllele> clippedReadCountsBuilder;
    private final ImmutableMultiset.Builder<PileAllele> lastPositionReadCountsBuilder;
    private final ImmutableMultiset.Builder<PileAllele> apparentMismapReadCountsBuilder;
    private final ImmutableMultiset.Builder<PileAllele> unmappedMateCountsBuilder;

    public Builder(GenomePosition position) {
      this.position = position;
      queriedAlleles = generateQueriedAlleles(position);
      basePilesBuilder = ImmutableSetMultimap.builder();
      weightedDepth = Maps.newHashMap();
      clippedReadCountsBuilder = ImmutableMultiset.builder();
      lastPositionReadCountsBuilder = ImmutableMultiset.builder();
      apparentMismapReadCountsBuilder = ImmutableMultiset.builder();
      unmappedMateCountsBuilder = ImmutableMultiset.builder();
    }

    public Builder addRecord(SAMRecord samRecord) {
      if (!samRecord.getDuplicateReadFlag()) {
        int readPos = samRecord.getReadPositionAtReferencePosition(position.getPosition()) - 1;
        if (readPos != -1) {
          PileAllele allele =
              queriedAlleles
                  .stream()
                  .filter(a -> a.supported(samRecord, readPos))
                  .findFirst()
                  .orElseGet(() -> getAppropriateAllele(samRecord, readPos));
          basePilesBuilder.put(allele, samRecord.hashCode());
          boolean countWeight = true;
          if (samRecord.getCigar().isClipped()) {
            clippedReadCountsBuilder.add(allele);
            countWeight = false;
          } else if (calcPercentReadMatchesRef(samRecord) >= MIN_PERCENT_BASES_MATCH) {
            apparentMismapReadCountsBuilder.add(allele);
            countWeight = false;
          }
          if (samRecord.getAlignmentStart() == position.getPosition()
              || samRecord.getAlignmentEnd() == position.getPosition()) {
            lastPositionReadCountsBuilder.add(allele);
          }
          if (samRecord.getMateUnmappedFlag()) {
            unmappedMateCountsBuilder.add(allele);
            countWeight = false;
          }
          if (countWeight) {
            weightedDepth.put(
                allele,
                weightedDepth.getOrDefault(allele, 0.0) + allele.weightedDepth(samRecord, readPos));
          }
        }
      }
      return this;
    }

    public Builder addAll(Stream<SAMRecord> samRecords) {
      samRecords.forEach(this::addRecord);
      return this;
    }

    public Pileup build() {
      return new Pileup(this);
    }

    private double calcPercentReadMatchesRef(SAMRecord samRecord) {
      return samRecord
              .getCigar()
              .getCigarElements()
              .stream()
              .filter(c -> c.getOperator().equals(CigarOperator.EQ))
              .mapToInt(CigarElement::getLength)
              .sum()
          / (double) samRecord.getReadLength();
    }
  }

  private static final double MIN_PERCENT_BASES_MATCH = 0.5;

  private final ImmutableSetMultimap<PileAllele, Integer> basePiles;
  private final ImmutableMap<PileAllele, Double> weightedBaseCounts;
  private final ImmutableMultiset<PileAllele> clippedReadCounts;
  private final ImmutableMultiset<PileAllele> lastPositionReadCounts;
  private final ImmutableMultiset<PileAllele> apparentMismapReadCounts;
  private final ImmutableMultiset<PileAllele> unmappedMateCounts;

  private final GenomePosition position;

  private Optional<Depth> depth = Optional.absent();

  public Pileup(Stream<SAMRecord> queriedRecords, GenomePosition position) {
    this(new Builder(position).addAll(queriedRecords));
  }

  private Pileup(Builder builder) {
    basePiles = builder.basePilesBuilder.build();
    weightedBaseCounts =
        ImmutableMap.<PileAllele, Double>builderWithExpectedSize(builder.weightedDepth.size())
            .putAll(builder.weightedDepth)
            .orderEntriesByValue(Comparator.reverseOrder())
            .build();
    clippedReadCounts = builder.clippedReadCountsBuilder.build();
    lastPositionReadCounts = builder.lastPositionReadCountsBuilder.build();
    apparentMismapReadCounts = builder.apparentMismapReadCountsBuilder.build();
    unmappedMateCounts = builder.unmappedMateCountsBuilder.build();
    position = builder.position;
  }

  private static PileAllele getAppropriateAllele(SAMRecord samRecord, int readPos) {
    byte base = samRecord.getReadBases()[readPos];
    return SNPAllele.of(base);
  }

  private static ImmutableList<PileAllele> generateQueriedAlleles(GenomePosition pos) {
    if (!(pos instanceof ReferencePosition)) return ImmutableList.of();
    ReferencePosition refPos = (ReferencePosition) pos;
    ImmutableList.Builder<PileAllele> queriedAllelesBuilder =
        ImmutableList.builderWithExpectedSize(2);
    queriedAllelesBuilder.add(refPos.getRefAllele());
    refPos.getAltAllele().toJavaUtil().ifPresent(queriedAllelesBuilder::add);
    return queriedAllelesBuilder.build();
  }

  /** @return Multiset of {@link PileAllele} counts */
  public ImmutableMultiset<PileAllele> getBaseCounts() {
    return basePiles.keys();
  }

  /**
   * @return Map from {@link PileAllele} to fraction of total depth for that {@link PileAllele},
   *     iteration order is in descending order of base fraction
   */
  public ImmutableMap<PileAllele, Double> getBaseFractions() {

    return ImmutableMap.copyOf(
        Maps.<PileAllele, Double>asMap(
            getBaseCounts().elementSet(),
            b -> getBaseCounts().count(b) / (double) getBaseCounts().size()));
  }

  /**
   * @return Map from {@link PileAllele} to weighted depth for that PileAllele, iteration order is
   *     in descending order of weighted base counts
   */
  public ImmutableMap<PileAllele, Double> getWeightedBaseCounts() {
    return weightedBaseCounts;
  }

  /** @return Multimap from {@link PileAllele} to index for the piled read */
  public ImmutableSetMultimap<PileAllele, Integer> getRecordsByBase() {
    return basePiles;
  }

  /** @return the clippedReadCounts */
  public ImmutableMultiset<PileAllele> getClippedReadCounts() {
    return clippedReadCounts;
  }

  /** @return the lastPositionReadCounts */
  public ImmutableMultiset<PileAllele> getLastPositionReadCounts() {
    return lastPositionReadCounts;
  }

  /** @return the apparentMismapReadCounts */
  public ImmutableMultiset<PileAllele> getApparentMismapReadCounts() {
    return apparentMismapReadCounts;
  }

  /** @return the unmappedMateCounts */
  public ImmutableMultiset<PileAllele> getUnmappedMateCounts() {
    return unmappedMateCounts;
  }

  /**
   * @return Map from {@link PileAllele} to weighted fraction of total weighted depth for that
   *     {@link PileAllele}, iteration order is in descending order of weighted base fraction
   */
  public ImmutableMap<PileAllele, Double> getWeightedBaseFractions() {
    final double weightedDepth =
        weightedBaseCounts.values().stream().mapToDouble(Double::valueOf).sum();
    return ImmutableMap.copyOf(Maps.transformValues(weightedBaseCounts, c -> c / weightedDepth));
  }

  /** @return the position */
  public GenomePosition getPosition() {
    return position;
  }

  private Depth setDepth() {
    depth = Optional.of(new Depth(this));
    return depth.get();
  }

  public Depth getDepth() {
    return depth.or(this::setDepth);
  }

  @Override
  public String toString() {
    return "Pileup [basePiles="
        + basePiles
        + ", weightedBaseCounts="
        + weightedBaseCounts
        + ", clippedReadCounts="
        + clippedReadCounts
        + ", lastPositionReadCounts="
        + lastPositionReadCounts
        + ", apparentMismapReadCounts="
        + apparentMismapReadCounts
        + ", unmappedMateCounts="
        + unmappedMateCounts
        + ", position="
        + position
        + ", depth="
        + depth
        + "]";
  }
}
