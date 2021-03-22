package org.pankratzlab.supernovo;

import java.io.Serializable;
import com.google.common.base.Optional;

public class GenomePosition implements Comparable<GenomePosition>, Serializable {

  private static final long serialVersionUID = 1L;

  protected final String contig;
  protected final int position;

  /**
   * @param contig
   * @param position
   */
  public GenomePosition(String contig, int position) {
    super();
    this.contig = contig;
    this.position = position;
  }

  /** @return the contig */
  public String getContig() {
    return contig;
  }

  /** @return the position */
  public int getPosition() {
    return position;
  }

  private Optional<Integer> locateChrNumber() {
    int i;
    for (i = 0; i < contig.length(); i++) {
      if (Character.isDigit(contig.charAt(i))) break;
    }
    StringBuilder digitBuilder = new StringBuilder();
    for (; i < contig.length(); i++) {
      char charAt = contig.charAt(i);
      if (Character.isDigit(charAt)) digitBuilder.append(charAt);
    }
    try {
      return Optional.of(Integer.parseInt(digitBuilder.toString()));
    } catch (NumberFormatException e) {
      return Optional.absent();
    }
  }

  @Override
  public int compareTo(GenomePosition o) {
    int contigCmp = compareContigs(o);
    if (contigCmp != 0) return contigCmp;
    int posCmp = Integer.compare(position, o.position);
    return posCmp;
  }

  private int compareContigs(GenomePosition o) {
    if (contig.equals(o.contig)) return 0;
    Optional<Integer> chrNumeric = locateChrNumber();
    Optional<Integer> otherChrNumeric = o.locateChrNumber();
    if (chrNumeric.isPresent() && otherChrNumeric.isPresent())
      return chrNumeric.get().compareTo(otherChrNumeric.get());
    if (chrNumeric.isPresent()) return -1;
    if (otherChrNumeric.isPresent()) return 1;
    return contig.compareTo(o.contig);
  }

  /* (non-Javadoc)
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((contig == null) ? 0 : contig.hashCode());
    result = prime * result + position;
    return result;
  }

  /* (non-Javadoc)
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj instanceof GenomePosition)) return false;
    GenomePosition other = (GenomePosition) obj;
    if (contig == null) {
      if (other.contig != null) return false;
    } else if (!contig.equals(other.contig)) return false;
    if (position != other.position) return false;
    return true;
  }

  @Override
  public String toString() {
    return "GenomePosition [contig=" + contig + ", position=" + position + "]";
  }
}
