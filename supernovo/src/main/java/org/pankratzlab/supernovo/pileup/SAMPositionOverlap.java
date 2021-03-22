package org.pankratzlab.supernovo.pileup;

import java.util.stream.Stream;
import htsjdk.samtools.SAMRecord;

public interface SAMPositionOverlap {

  /** @return the records */
  Stream<SAMRecord> getRecords();
}
