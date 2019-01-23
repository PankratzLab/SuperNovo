package org.pankratzlab.supernovo.pileup;

import java.util.Comparator;
import java.util.Map;
import org.pankratzlab.supernovo.Position;
import org.pankratzlab.supernovo.utilities.Phred;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.Maps;
import htsjdk.samtools.SAMRecord;

public class SAMRecordPileup extends AbstractPileup {

  private final ImmutableSetMultimap<Byte, Integer> basePiles;
  private final ImmutableMap<Byte, Double> weightedBaseCounts;
  private final ImmutableList<SAMRecord> queriedRecords;
  private final ImmutableMultiset<Byte> clippedReadCounts;
  private final ImmutableMultiset<Byte> unmappedMateCounts;

  public SAMRecordPileup(ImmutableList<SAMRecord> queriedRecords, Position position) {
    super();
    ImmutableSetMultimap.Builder<Byte, Integer> basePilesBuilder = ImmutableSetMultimap.builder();
    Map<Byte, Double> weightedDepth = Maps.newHashMap();
    ImmutableMultiset.Builder<Byte> clippedReadCountsBuilder = ImmutableMultiset.builder();
    ImmutableMultiset.Builder<Byte> unmappedMateCountsBuilder = ImmutableMultiset.builder();
    for (int i = 0; i < queriedRecords.size(); i++) {
      SAMRecord samRecord = queriedRecords.get(i);
      int readPos = samRecord.getReadPositionAtReferencePosition(position.getPosition()) - 1;
      if (readPos != -1) {
        Byte base = samRecord.getReadBases()[readPos];
        basePilesBuilder.put(base, i);
        double accuracy =
            Phred.getAccuracy(samRecord.getBaseQualities()[readPos])
                * Phred.getAccuracy(samRecord.getMappingQuality());
        weightedDepth.put(base, weightedDepth.getOrDefault(base, 0.0) + accuracy);
        if (samRecord.getCigar().isClipped()) clippedReadCountsBuilder.add(base);
        if (samRecord.getMateUnmappedFlag()) unmappedMateCountsBuilder.add(base);
      }
    }
    basePiles = basePilesBuilder.build();
    weightedBaseCounts =
        ImmutableMap.<Byte, Double>builderWithExpectedSize(weightedDepth.size())
            .putAll(weightedDepth)
            .orderEntriesByValue(Comparator.reverseOrder())
            .build();
    clippedReadCounts = clippedReadCountsBuilder.build();
    unmappedMateCounts = unmappedMateCountsBuilder.build();
    this.queriedRecords = queriedRecords;
  }

  @Override
  public ImmutableMultiset<Byte> getBaseCounts() {
    return basePiles.keys();
  }

  @Override
  public ImmutableMap<Byte, Double> getWeightedBaseCounts() {
    return weightedBaseCounts;
  }

  @Override
  public ImmutableSetMultimap<Byte, Integer> getRecordsByBase() {
    return basePiles;
  }

  @Override
  public ImmutableList<SAMRecord> getRecords() {
    return queriedRecords;
  }

  /** @return the clippedReadCounts */
  @Override
  public ImmutableMultiset<Byte> getClippedReadCounts() {
    return clippedReadCounts;
  }

  /** @return the unmappedMateCounts */
  @Override
  public ImmutableMultiset<Byte> getUnmappedMateCounts() {
    return unmappedMateCounts;
  }

  @Override
  public ImmutableMap<Byte, Double> getWeightedBaseFractions() {
    final double weightedDepth =
        weightedBaseCounts.values().stream().mapToDouble(Double::valueOf).sum();
    return ImmutableMap.copyOf(Maps.transformValues(weightedBaseCounts, c -> c / weightedDepth));
  }
}
