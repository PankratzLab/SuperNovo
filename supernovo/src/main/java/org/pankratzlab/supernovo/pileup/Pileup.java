package org.pankratzlab.supernovo.pileup;

import java.io.Serializable;
import java.util.Comparator;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;
import org.pankratzlab.supernovo.GenomePosition;
import org.pankratzlab.supernovo.PileAllele;
import org.pankratzlab.supernovo.ReferencePosition;
import org.pankratzlab.supernovo.SNPAllele;
import com.google.common.base.Optional;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimaps;
import com.google.common.collect.Sets;
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

  /** @return Multimap from {@link PileAllele} to hash code for the piled read */
  public ImmutableSetMultimap<PileAllele, Integer> getRecordsByBase() {
    return basePiles;
  }

  /** @return Set of hash code for each piled read */
  public Set<Integer> getRecords() {
    return Multimaps.asMap(basePiles)
        .values()
        .stream()
        .reduce(Sets::union)
        .orElse(ImmutableSet.of());
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

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result =
        prime * result
            + ((apparentMismapReadCounts == null) ? 0 : apparentMismapReadCounts.hashCode());
    result = prime * result + ((basePiles == null) ? 0 : basePiles.hashCode());
    result = prime * result + ((clippedReadCounts == null) ? 0 : clippedReadCounts.hashCode());
    result = prime * result + ((depth == null) ? 0 : depth.hashCode());
    result =
        prime * result + ((lastPositionReadCounts == null) ? 0 : lastPositionReadCounts.hashCode());
    result = prime * result + ((position == null) ? 0 : position.hashCode());
    result = prime * result + ((unmappedMateCounts == null) ? 0 : unmappedMateCounts.hashCode());
    result = prime * result + ((weightedBaseCounts == null) ? 0 : weightedBaseCounts.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj instanceof Pileup)) return false;
    Pileup other = (Pileup) obj;
    if (apparentMismapReadCounts == null) {
      if (other.apparentMismapReadCounts != null) return false;
    } else if (!apparentMismapReadCounts.equals(other.apparentMismapReadCounts)) return false;
    if (basePiles == null) {
      if (other.basePiles != null) return false;
    } else if (!basePiles.equals(other.basePiles)) return false;
    if (clippedReadCounts == null) {
      if (other.clippedReadCounts != null) return false;
    } else if (!clippedReadCounts.equals(other.clippedReadCounts)) return false;
    if (depth == null) {
      if (other.depth != null) return false;
    } else if (!depth.equals(other.depth)) return false;
    if (lastPositionReadCounts == null) {
      if (other.lastPositionReadCounts != null) return false;
    } else if (!lastPositionReadCounts.equals(other.lastPositionReadCounts)) return false;
    if (position == null) {
      if (other.position != null) return false;
    } else if (!position.equals(other.position)) return false;
    if (unmappedMateCounts == null) {
      if (other.unmappedMateCounts != null) return false;
    } else if (!unmappedMateCounts.equals(other.unmappedMateCounts)) return false;
    if (weightedBaseCounts == null) {
      if (other.weightedBaseCounts != null) return false;
    } else if (!weightedBaseCounts.equals(other.weightedBaseCounts)) return false;
    return true;
  }
}
