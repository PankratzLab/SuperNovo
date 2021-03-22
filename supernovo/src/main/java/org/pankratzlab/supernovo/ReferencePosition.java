package org.pankratzlab.supernovo;

import java.io.Serializable;
import com.google.common.base.Optional;
import com.google.common.collect.ImmutableList;
import com.google.common.primitives.Bytes;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class ReferencePosition extends GenomePosition implements Serializable {

  /** */
  private static final long serialVersionUID = 1L;

  private final PileAllele refAllele;
  private final Optional<PileAllele> altAllele;

  /**
   * @param contig
   * @param position
   * @param refAllele
   * @param altAllele
   */
  private ReferencePosition(
      String contig, int position, PileAllele refAllele, Optional<PileAllele> altAllele) {
    super(contig, position);
    this.refAllele = refAllele;
    this.altAllele = altAllele;
  }

  /**
   * @param contig
   * @param position
   * @param refAllele
   * @param altAllele
   */
  public ReferencePosition(
      String contig, int position, PileAllele refAllele, PileAllele altAllele) {
    this(contig, position, refAllele, Optional.of(altAllele));
  }

  /**
   * @param contig
   * @param position
   * @param refAllele
   * @param altAllele
   */
  public ReferencePosition(String contig, int position, PileAllele refAllele) {
    this(contig, position, refAllele, Optional.absent());
  }

  public static ReferencePosition fromVariantContext(VariantContext vc, Allele ref, Allele alt) {
    final PileAllele refAllele;
    final PileAllele altAllele;
    if (ref.length() == 1 && alt.length() == 1) {
      refAllele = SNPAllele.of(ref.getBases()[0]);
      altAllele = SNPAllele.of(alt.getBases()[0]);
    } else if (ref.length() == 1 && alt.length() > 1) {
      altAllele = generateInsertionAllele(alt, ref);
      refAllele = ((InsertionAllele) altAllele).getNonInsertionAllele();
    } else if (alt.length() == 1 && ref.length() > 1) {
      refAllele = generateInsertionAllele(ref, alt);
      altAllele = ((InsertionAllele) refAllele).getNonInsertionAllele();
    } else throw new IllegalArgumentException("Only SNPs and Indels are supported");
    return new ReferencePosition(vc.getContig(), vc.getStart(), refAllele, altAllele);
  }

  private static InsertionAllele generateInsertionAllele(Allele ins, Allele del) {
    byte preBase = del.getBases()[0];
    if (preBase != ins.getBases()[0])
      throw new IllegalArgumentException(
          "Indels must match on first base (ins: "
              + ins.toString()
              + ", del: "
              + del.toString()
              + ")");
    return new InsertionAllele(
        SNPAllele.of(preBase),
        ImmutableList.copyOf(Bytes.asList(ins.getBases()).subList(1, ins.getBases().length)));
  }

  /** @return the refAllele */
  public PileAllele getRefAllele() {
    return refAllele;
  }

  /** @return the altAllele */
  public Optional<PileAllele> getAltAllele() {
    return altAllele;
  }

  /* (non-Javadoc)
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((altAllele == null) ? 0 : altAllele.hashCode());
    result = prime * result + ((contig == null) ? 0 : contig.hashCode());
    result = prime * result + position;
    result = prime * result + ((refAllele == null) ? 0 : refAllele.hashCode());
    return result;
  }

  /* (non-Javadoc)
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj instanceof ReferencePosition)) return false;
    ReferencePosition other = (ReferencePosition) obj;
    if (altAllele == null) {
      if (other.altAllele != null) return false;
    } else if (!altAllele.equals(other.altAllele)) return false;
    if (contig == null) {
      if (other.contig != null) return false;
    } else if (!contig.equals(other.contig)) return false;
    if (position != other.position) return false;
    if (refAllele == null) {
      if (other.refAllele != null) return false;
    } else if (!refAllele.equals(other.refAllele)) return false;
    return true;
  }
}
