package org.pankratzlab.supernovo;

import java.io.Serializable;
import htsjdk.samtools.SAMRecord;

public interface PileAllele extends Serializable {

  /**
   * @param record to test
   * @param readPos to query
   * @return true if this read supports this {@link PileAllele}
   */
  boolean supported(SAMRecord record, int readPos);

  /**
   * @param record to test
   * @param readPos to query
   * @return true if support for this {@link PileAllele} is based on a clipped portion of record
   */
  boolean clipped(SAMRecord record, int readPos);

  /**
   * @param record to test
   * @param readPos to query
   * @return weighted depth for allele
   */
  double weightedDepth(SAMRecord record, int readPos);

  /** @return String representation of the allele */
  @Override
  String toString();

  @Override
  boolean equals(Object obj);

  @Override
  int hashCode();
}
