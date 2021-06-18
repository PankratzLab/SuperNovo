package org.pankratzlab.supernovo;

import java.io.File;
import org.pankratzlab.supernovo.output.DeNovoResult;
import org.pankratzlab.supernovo.pileup.Pileup;
import org.pankratzlab.supernovo.pileup.PileupCache;
import com.google.common.collect.ImmutableList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class TrioEvaluator extends AbstractEvaluator {

  private final PileupCache p1Pileups;
  private final PileupCache p2Pileups;

  public TrioEvaluator(File childBam, File p1Bam, File p2Bam) {
    super(childBam);
    this.p1Pileups = new PileupCache(p1Bam);
    this.p2Pileups = new PileupCache(p2Bam);
  }

  @Override
  protected boolean keepVariant(VariantContext vc) {
    return super.keepVariant(vc) && !seenInParentVCF(vc);
  }

  @Override
  protected DeNovoResult generateDeNovoResult(
      ReferencePosition pos, Pileup childPile, PileupCache childPileups) {
    return new DeNovoResult(
        pos,
        new HaplotypeEvaluator(childPile, childPileups, p1Pileups, p2Pileups, this)
            .haplotypeConcordance(),
        generateSample(App.getInstance().getChildID(), pos, childPile, childPile),
        generateSample(App.getInstance().getP1ID(), pos, p1Pileups.get(pos), childPile),
        generateSample(App.getInstance().getP2ID(), pos, p2Pileups.get(pos), childPile));
  }

  private boolean seenInParentVCF(VariantContext vc) {
    Allele altAllele = getPossibleDnAlt(vc);
    for (String parentID :
        ImmutableList.of(App.getInstance().getP1ID(), App.getInstance().getP2ID())) {
      Genotype geno = vc.getGenotype(parentID);
      int parentAlleleIndex = geno.getAlleles().indexOf(altAllele);
      if (parentAlleleIndex >= 0
          && geno.getAD()[parentAlleleIndex] > App.getInstance().getVcfMaxParentAD()) {
        return true;
      }
    }
    return false;
  }
}
