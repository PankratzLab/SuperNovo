package org.pankratzlab.supernovo.pileup;

import java.util.stream.Stream;
import org.pankratzlab.supernovo.GenomePosition;
import htsjdk.samtools.SAMRecord;

public class SAMPositionRecordsOverlap implements SAMPositionOverlap {

  private final Stream<SAMRecord> records;

  public SAMPositionRecordsOverlap(Stream<SAMRecord> records, GenomePosition position) {
    this.records = records;
  }

  /* (non-Javadoc)
   * @see org.pankratzlab.supernovo.pileup.SAMPositionOverlap#getRecords()
   */
  @Override
  public Stream<SAMRecord> getRecords() {
    return records;
  }
}
