package org.pankratzlab.supernovo.pileup;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.stream.Stream;
import org.pankratzlab.supernovo.GenomePosition;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class SAMPositionQueryOverlap implements SAMPositionOverlap, Closeable {

  public static final SamReaderFactory SR_FACTORY = SamReaderFactory.make();

  private final Stream<SAMRecord> records;
  private final SAMRecordIterator iterator;
  private final SamReader reader;

  public SAMPositionQueryOverlap(File bam, GenomePosition position) {
    reader = SR_FACTORY.open(bam);
    iterator =
        reader.queryOverlapping(
            position.getContig(), position.getPosition(), position.getPosition());
    records = iterator.stream();
  }

  /* (non-Javadoc)
   * @see org.pankratzlab.supernovo.pileup.SAMPositionOverlap#getRecords()
   */
  @Override
  public Stream<SAMRecord> getRecords() {
    return records;
  }

  @Override
  public void close() {
    iterator.close();
    try {
      reader.close();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
}
