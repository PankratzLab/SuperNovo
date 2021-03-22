package org.pankratzlab.supernovo;

import java.util.stream.IntStream;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;

public class InsertionAllele extends AbstractPileAllele {

  /** */
  private static final long serialVersionUID = 1L;

  private class NonInsertionAllele extends AbstractPileAllele {

    /** */
    private static final long serialVersionUID = 1L;

    private NonInsertionAllele() {
      super(preInsertionBase.toString());
    }

    @Override
    public boolean supported(SAMRecord record, int readPos) {
      return InsertionAllele.this.supportType(record, readPos).equals(Support.NO_INSERTION);
    }

    @Override
    public boolean clipped(SAMRecord record, int readPos) {
      return false;
    }

    public InsertionAllele getInsertionAllele() {
      return InsertionAllele.this;
    }

    @Override
    public double weightedDepth(SAMRecord samRecord, int readPos) {
      return singlePosWeightedDepth(samRecord, readPos);
    }

    /* (non-Javadoc)
     * @see java.lang.Object#hashCode()
     */
    @Override
    public int hashCode() {
      final int prime = 31;
      int result = super.hashCode();
      result = prime * result + InsertionAllele.this.hashCode();
      return result;
    }

    /* (non-Javadoc)
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (!super.equals(obj)) return false;
      if (!(obj instanceof NonInsertionAllele)) return false;
      NonInsertionAllele other = (NonInsertionAllele) obj;
      return getInsertionAllele().equals(other.getInsertionAllele());
    }
  }

  private enum Support {
    CLIPPED_INSERTION(true),
    CLEAN_INSERTION(true),
    NO_INSERTION(false),
    OTHER(false);

    private final boolean supported;

    /** @param supported */
    private Support(boolean supported) {
      this.supported = supported;
    }
  }

  private final SNPAllele preInsertionBase;
  private final ImmutableList<Byte> insertedBases;
  private final NonInsertionAllele nonInsertionAllele;

  /**
   * @param preInsertionBase
   * @param insertedBases
   */
  public InsertionAllele(SNPAllele preInsertionBase, ImmutableList<Byte> insertedBases) {
    super(alleleString(preInsertionBase, insertedBases));
    this.preInsertionBase = preInsertionBase;
    this.insertedBases = insertedBases;
    this.nonInsertionAllele = new NonInsertionAllele();
  }

  @Override
  public boolean supported(final SAMRecord record, final int readPos) {
    return supportType(record, readPos).supported;
  }

  @Override
  public boolean clipped(SAMRecord record, int readPos) {
    return supportType(record, readPos).equals(Support.CLIPPED_INSERTION);
  }

  @Override
  public double weightedDepth(SAMRecord samRecord, int readPos) {
    int length = insertedBases.size() + 1;
    byte[] readBases = samRecord.getReadBases();
    int limit = Integer.min(readBases.length, readPos + length);
    return IntStream.range(readPos, limit)
        .mapToDouble(i -> singlePosWeightedDepth(samRecord, i))
        .average()
        .orElseGet(() -> Double.valueOf(0.0));
  }

  private Support supportType(final SAMRecord record, final int readPos) {
    if (preInsertionBase.supported(record, readPos)) {
      boolean support = false;
      boolean clipped = false;
      byte[] readBases = record.getReadBases();
      int offset = readPos + 1;
      for (int i = 0; i < insertedBases.size() && i + offset < readBases.length; i++) {
        if (readBases[i + offset] == insertedBases.get(i)) {
          support = true;
          if (record.getReferencePositionAtReadPosition(i + offset) == 0) clipped = true;
        } else {
          if (support) return Support.OTHER;
          return Support.NO_INSERTION;
        }
      }
      if (support) return clipped ? Support.CLIPPED_INSERTION : Support.CLEAN_INSERTION;
      return Support.OTHER;
    }
    return Support.OTHER;
  }

  private static String alleleString(
      SNPAllele preInsertionBase, ImmutableList<Byte> insertedBases) {
    String preInsertion = preInsertionBase.toString();
    StringBuilder alleleBuilder = new StringBuilder(preInsertion.length() + insertedBases.size());
    alleleBuilder.append(preInsertion);
    insertedBases.forEach(b -> alleleBuilder.append((char) b.byteValue()));
    return alleleBuilder.toString();
  }

  /** @return a {@link PileAllele} that represents the alternative to this insertion */
  public NonInsertionAllele getNonInsertionAllele() {
    return nonInsertionAllele;
  }

  /* (non-Javadoc)
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode() {
    final int prime = 31;
    int result = super.hashCode();
    result = prime * result + ((insertedBases == null) ? 0 : insertedBases.hashCode());
    result = prime * result + ((preInsertionBase == null) ? 0 : preInsertionBase.hashCode());
    return result;
  }

  /* (non-Javadoc)
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (!super.equals(obj)) return false;
    if (!(obj instanceof InsertionAllele)) return false;
    InsertionAllele other = (InsertionAllele) obj;
    if (insertedBases == null) {
      if (other.insertedBases != null) return false;
    } else if (!insertedBases.equals(other.insertedBases)) return false;
    if (preInsertionBase == null) {
      if (other.preInsertionBase != null) return false;
    } else if (!preInsertionBase.equals(other.preInsertionBase)) return false;
    return true;
  }
}
