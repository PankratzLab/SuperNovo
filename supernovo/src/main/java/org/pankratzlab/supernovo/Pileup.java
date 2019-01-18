package org.pankratzlab.supernovo;

import org.pankratzlab.supernovo.SAMRecordPileup.PiledRecord;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.ImmutableMultiset;

public interface Pileup {

  /**
   * @return Map from byte value of base to weighted depth for that base, iteration order is in
   *     descending order of weighted base counts
   */
  public ImmutableMap<Byte, Double> getWeightedBaseCounts();

  /** @return Multiset of byte values of bases */
  public ImmutableMultiset<Byte> getBaseCounts();

  /** @return Multimap from byte value of bases to {@link PiledRecord}s for that base */
  public ImmutableMultimap<Byte, PiledRecord> getRecordsByBase();
}
