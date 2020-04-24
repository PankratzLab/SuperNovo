package org.pankratzlab.supernovo;

import java.io.Serializable;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.Function;
import org.pankratzlab.supernovo.output.DeNovoResult;
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

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + adjacentDeNovos;
      result = prime * result + ((concordances == null) ? 0 : concordances.hashCode());
      result = prime * result + otherBiallelics;
      result = prime * result + otherDeNovos;
      result = prime * result + otherTriallelics;
      result = prime * result + otherVariants;
      return result;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (!(obj instanceof Result)) return false;
      Result other = (Result) obj;
      if (adjacentDeNovos != other.adjacentDeNovos) return false;
      if (concordances == null) {
        if (other.concordances != null) return false;
      } else if (!concordances.equals(other.concordances)) return false;
      if (otherBiallelics != other.otherBiallelics) return false;
      if (otherDeNovos != other.otherDeNovos) return false;
      if (otherTriallelics != other.otherTriallelics) return false;
      if (otherVariants != other.otherVariants) return false;
      return true;
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
        if (TrioEvaluator.looksVariant(searchPileup.getDepth())) {
          otherVariants++;
          if (TrioEvaluator.moreThanTwoViableAlleles(searchPileup)) {
            otherTriallelics++;
          } else {
            otherBiallelics++;
            concordance(childPile, searchPileup).ifPresent(concordances::add);
          }
        }
        if (TrioEvaluator.passesAllelicFrac(searchPileup.getDepth())
            && concordance(childPile, searchPileup).orElse(0.0)
                >= DeNovoResult.MIN_HAPLOTYPE_CONCORDANCE
            && TrioEvaluator.looksDenovo(
                searchPileup, p1Piles.apply(searchPosition), p2Piles.apply(searchPosition))) {
          otherDenovoPositions.add(searchPos);
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

  private static OptionalDouble concordance(Pileup base, Pileup search) {
    Set<Integer> h1 = base.getDepth().allelicRecords(Allele.A1);
    Set<Integer> h2 = base.getDepth().allelicRecords(Allele.A2);

    Set<Integer> search1 = search.getDepth().allelicRecords(Allele.A1);
    Set<Integer> search2 = search.getDepth().allelicRecords(Allele.A2);

    Set<Integer> h1Overlap = Sets.intersection(h1, search.getRecords());
    Set<Integer> h2Overlap = Sets.intersection(h2, search.getRecords());

    int h1Size = h1Overlap.size();
    int h2Size = h2Overlap.size();

    if (h1Size == 0 && h2Size == 0) return OptionalDouble.empty();

    final double h1CisConc =
        h1Size == 0 ? 1.0 : Sets.intersection(h1, search1).size() / (double) h1Size;
    final double h2CisConc =
        h2Size == 0 ? 1.0 : Sets.intersection(h2, search2).size() / (double) h2Size;

    double cisConc = Math.min(h1CisConc, h2CisConc);

    final double h1TransConc =
        h1Size == 0 ? 1.0 : Sets.intersection(h1, search2).size() / (double) h1Size;
    final double h2TransConc =
        h2Size == 0 ? 1.0 : Sets.intersection(h2, search1).size() / (double) h2Size;

    double transConc = Math.min(h1TransConc, h2TransConc);

    return OptionalDouble.of(Math.max(cisConc, transConc));
  }
}
