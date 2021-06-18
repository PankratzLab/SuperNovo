package org.pankratzlab.supernovo;

import java.io.File;
import org.pankratzlab.supernovo.output.DeNovoResult;
import org.pankratzlab.supernovo.pileup.Pileup;
import org.pankratzlab.supernovo.pileup.PileupCache;

public class VariantEvaluator extends AbstractEvaluator {

  public VariantEvaluator(File bam) {
    super(bam);
  }

  @Override
  protected DeNovoResult generateDeNovoResult(
      ReferencePosition pos, Pileup childPile, PileupCache childPileups) {
    return new DeNovoResult(
        pos,
        new HaplotypeEvaluator(childPile, childPileups, this).haplotypeConcordance(),
        generateSample(App.getInstance().getChildID(), pos, childPile, childPile));
  }
}
