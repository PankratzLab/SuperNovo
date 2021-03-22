package org.pankratzlab.supernovo;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import htsjdk.samtools.SAMRecord;

public class SNPAllele extends AbstractPileAllele {

  /** */
  private static final long serialVersionUID = 1L;

  private static LoadingCache<Byte, SNPAllele> cache =
      CacheBuilder.newBuilder().softValues().build(CacheLoader.from(SNPAllele::new));

  public static final SNPAllele A = of((byte) 'A');
  public static final SNPAllele T = of((byte) 'T');
  public static final SNPAllele C = of((byte) 'C');
  public static final SNPAllele G = of((byte) 'G');

  private final byte base;

  /** @param base */
  private SNPAllele(byte base) {
    super(String.valueOf((char) base));
    this.base = base;
  }

  public static SNPAllele of(byte base) {
    return cache.getUnchecked(base);
  }

  @Override
  public boolean supported(SAMRecord record, int readPos) {
    return readPos != -1 && record.getReadBases()[readPos] == base;
  }

  @Override
  public boolean clipped(SAMRecord record, int readPos) {
    return record.getReferencePositionAtReadPosition(readPos) != 0;
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
    result = prime * result + base;
    return result;
  }

  /* (non-Javadoc)
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (!super.equals(obj)) return false;
    if (!(obj instanceof SNPAllele)) return false;
    SNPAllele other = (SNPAllele) obj;
    if (base != other.base) return false;
    return true;
  }
}
