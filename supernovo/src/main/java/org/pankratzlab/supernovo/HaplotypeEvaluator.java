package org.pankratzlab.supernovo;

import java.io.Serializable;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import org.pankratzlab.supernovo.pileup.Depth.Allele;
import org.pankratzlab.supernovo.pileup.Pileup;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Sets;

public class HaplotypeEvaluator {

  public static class Result implements Serializable {
    /** */
    private static final long serialVersionUID = 1L;

    private final int otherVariants;
    private final int otherTriallelics;
    private final int otherBiallelics;
    private final int adjacentDeNovos;
    private final int otherDeNovos;
    private final List<Double> concordances;
    /**
     * @param otherVariants
     * @param otherTriallelics
     * @param otherBiallelics
     * @param adjacentDeNovos TODO
     * @param otherDeNovos
     * @param concordances
     * @param adjacentDeNovos
     */
    public Result(
        int otherVariants,
        int otherTriallelics,
        int otherBiallelics,
        int adjacentDeNovos,
        int otherDeNovos,
        List<Double> concordances) {
      super();
      this.otherVariants = otherVariants;
      this.otherTriallelics = otherTriallelics;
      this.otherBiallelics = otherBiallelics;
      this.adjacentDeNovos = adjacentDeNovos;
      this.otherDeNovos = otherDeNovos;
      this.concordances = concordances;
    }
    /** @return the otherVariants */
    public int getOtherVariants() {
      return otherVariants;
    }
    /** @return the otherTriallelics */
    public int getOtherTriallelics() {
      return otherTriallelics;
    }
    /** @return the otherBiallelics */
    public int getOtherBiallelics() {
      return otherBiallelics;
    }
    /** @return the adjacentDeNovos */
    public int getAdjacentDeNovos() {
      return adjacentDeNovos;
    }
    /** @return the otherDeNovos */
    public int getOtherDeNovos() {
      return otherDeNovos;
    }
    /** @return the concordances */
    public List<Double> getConcordances() {
      return concordances;
    }
  }

  private static final int HAPLOTYPE_SEARCH_DISTANCE = 150;

  private final Pileup childPile;
  private final GenomePosition pos;
  private final Function<GenomePosition, Pileup> childPiles;
  private final Function<GenomePosition, Pileup> p1Piles;
  private final Function<GenomePosition, Pileup> p2Piles;
  /**
   * @param child
   * @param p1
   * @param p2
   */
  public HaplotypeEvaluator(
      Pileup childPile,
      Function<GenomePosition, Pileup> childPiles,
      Function<GenomePosition, Pileup> p1Piles,
      Function<GenomePosition, Pileup> p2Piles) {
    super();
    this.childPile = childPile;
    this.pos = childPile.getPosition();
    this.childPiles = childPiles;
    this.p1Piles = p1Piles;
    this.p2Piles = p2Piles;
  }

  public Result haplotypeConcordance() {
    int startSearch = Integer.max(0, pos.getPosition() - HAPLOTYPE_SEARCH_DISTANCE);
    int stopSearch = pos.getPosition() + HAPLOTYPE_SEARCH_DISTANCE;

    Set<Integer> otherDenovoPositions = Sets.newHashSet();
    int otherTriallelics = 0;
    int otherBiallelics = 0;
    int otherVariants = 0;
    ImmutableList.Builder<Double> concordances = ImmutableList.builder();

    for (int searchPos = startSearch; searchPos < stopSearch; searchPos++) {
      if (searchPos == pos.getPosition()) continue;
      GenomePosition searchPosition = new GenomePosition(pos.getContig(), searchPos);
      Pileup searchPileup = childPiles.apply(searchPosition);
      if (searchPileup.getDepth().getBiAlleles().size() == 2) {
        otherVariants++;
        if (TrioEvaluator.moreThanTwoViableAlleles(searchPileup)) {
          otherTriallelics++;
        } else {
          otherBiallelics++;
          concordances.add(concordance(childPile, searchPileup));
          if (TrioEvaluator.looksDenovo(
              searchPileup, p1Piles.apply(searchPosition), p2Piles.apply(searchPosition))) {
            otherDenovoPositions.add(searchPos);
          }
        }
      }
    }
    int adjacentDeNovos = calculateAdjacentDenovos(otherDenovoPositions);
    int otherDenovos = otherDenovoPositions.size() - adjacentDeNovos;
    return new Result(
        otherVariants,
        otherTriallelics,
        otherBiallelics,
        adjacentDeNovos,
        otherDenovos,
        concordances.build());
  }

  private int calculateAdjacentDenovos(Set<Integer> otherDenovoPositions) {
    int adjacentDeNovos = 0;
    int searchPos = pos.getPosition();
    while (otherDenovoPositions.contains(++searchPos)) adjacentDeNovos++;
    searchPos = pos.getPosition();
    while (otherDenovoPositions.contains(--searchPos)) adjacentDeNovos++;
    return adjacentDeNovos;
  }

  private static double concordance(Pileup base, Pileup search) {
    Set<Integer> h1 = base.getDepth().allelicRecords(Allele.A1);
    Set<Integer> h2 = base.getDepth().allelicRecords(Allele.A2);

    Set<Integer> search1 = search.getDepth().allelicRecords(Allele.A1);
    Set<Integer> search2 = search.getDepth().allelicRecords(Allele.A2);

    double totalOverlap = search.getDepth().rawTotalDepth();
    int maxOverlap =
        Integer.max(
            Sets.intersection(search1, h1).size() + Sets.intersection(search2, h2).size(),
            Sets.intersection(search1, h2).size() + Sets.intersection(search2, h1).size());
    return maxOverlap / totalOverlap;
  }
}
