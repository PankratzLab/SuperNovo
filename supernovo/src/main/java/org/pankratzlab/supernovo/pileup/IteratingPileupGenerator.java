package org.pankratzlab.supernovo.pileup;

public class IteratingPileupGenerator {

  //  private final TreeMultimap<GenomePosition, SAMRecord> mappedRecords;
  //  private final Iterator<SAMRecord> records;
  //  private GenomePosition last;
  //
  //  public IteratingPileupGenerator(Iterator<SAMRecord> records) {
  //    this.records = records;
  //    mappedRecords = TreeMultimap.create(Ordering.natural(), new SAMRecordQueryNameComparator());
  //  }
  //
  //  public Pileup getPileup(GenomePosition position) {
  //    if (position.compareTo(last) > 0) {
  //      mappedRecords.asMap().headMap(toKey)
  //      last = position;
  //      while (records.hasNext()) {
  //        SAMRecord record = records.next();
  //      }
  //    }
  //  }
  //
  //  private void addRecord(SAMRecord record) {
  //    IntStream.range(record.getAlignmentStart(), record.getAlignmentEnd() + 1)
  //        .mapToObj(pos -> new GenomePosition(record.getContig(), pos))
  //        .forEach(gp -> mappedRecords.put(gp, record));
  //  }
}
