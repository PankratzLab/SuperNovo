package org.pankratzlab.supernovo.contam;

import java.util.SortedSet;
import java.util.stream.IntStream;
import com.google.common.collect.ConcurrentHashMultiset;
import com.google.common.collect.ImmutableSortedSet;
import com.google.common.collect.MoreCollectors;
import com.google.common.collect.Multiset;
import com.google.common.collect.Multisets;
import com.google.common.collect.Ordering;
import htsjdk.variant.variantcontext.Genotype;

public class AlleleRatio {

  private final SortedSet<Double> bins;
  private final Multiset<Double> altFracBinCounts;
  private final int minAltDepth;
  private final int minTotalDepth;

  public AlleleRatio(
      final ImmutableSortedSet<Double> bins, final int minAltDepth, final int minDepth) {
    this.bins = bins;
    this.minAltDepth = minAltDepth;
    this.minTotalDepth = minDepth;
    altFracBinCounts = ConcurrentHashMultiset.create();
  }

  public void addVariant(Genotype geno) {
    if (geno.hasAD()) {
      int[] ad = geno.getAD();
      int refCount = ad[0];
      int altCount = IntStream.range(1, ad.length).map(i -> ad[i]).max().orElse(0);
      int biallelicDepth = refCount + altCount;
      if (altCount >= minAltDepth && biallelicDepth >= minTotalDepth) {
        double altRatio = altCount / (double) biallelicDepth;
        altFracBinCounts.add(
            bins.tailSet(altRatio)
                .stream()
                .limit(1)
                .collect(MoreCollectors.toOptional())
                .orElseGet(bins::last));
      }
    }
  }

  /** @return the bins */
  public SortedSet<Double> getBins() {
    return bins;
  }

  /** @return the altFracBins */
  public Multiset<Double> getAltFracBinCounts() {
    return Multisets.unmodifiableMultiset(altFracBinCounts);
  }

  public static ImmutableSortedSet<Double> calculateBins(int nBins) {
    double binSize = 1.0 / nBins;
    return IntStream.range(1, nBins + 1)
        .mapToDouble(bin -> binSize * bin)
        .boxed()
        .collect(ImmutableSortedSet.toImmutableSortedSet(Ordering.natural()));
  }
}
