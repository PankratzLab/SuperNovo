package org.pankratzlab.supernovo.output;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.lang.ProcessBuilder.Redirect;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;
import org.apache.logging.log4j.Level;
import org.pankratzlab.supernovo.App;
import org.pankratzlab.supernovo.HaplotypeEvaluator;
import org.pankratzlab.supernovo.PileAllele;
import org.pankratzlab.supernovo.ReferencePosition;
import org.pankratzlab.supernovo.SNPAllele;
import org.pankratzlab.supernovo.TrioEvaluator;
import org.pankratzlab.supernovo.pileup.Depth;
import org.pankratzlab.supernovo.pileup.Pileup;
import com.google.common.base.Optional;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class DeNovoResult implements OutputFields, Serializable {

  /** */
  private static final long serialVersionUID = 1L;

  public static class Sample implements OutputFields, Serializable {

    /** */
    private static final long serialVersionUID = 1L;

    public final String id;
    public final int rawDepth;
    public final int refRawDepth;
    public final Optional<Integer> altRawDepth;
    public final int a1RawDepth;
    public final int a2RawDepth;
    public final int a_rawDepth;
    public final int t_rawDepth;
    public final int c_rawDepth;
    public final int g_rawDepth;
    public final int a1ClippedReads;
    public final int a2ClippedReads;
    public final int a1ApparentMismapReads;
    public final int a2ApparentMismapReads;
    public final int a1UnmappedMateReads;
    public final int a2UnmappedMateReads;
    public final double weightedDepth;
    public final double refWeightedDepth;
    public final Optional<Double> altWeightedDepth;
    public final double a1WeightedDepth;
    public final double a2WeightedDepth;
    public final double a_weightedDepth;
    public final double t_weightedDepth;
    public final double c_weightedDepth;
    public final double g_weightedDepth;

    private final Pileup pileup;

    /**
     * @param id
     * @param pileup
     * @param pos TODO
     */
    public Sample(
        String id,
        Pileup pileup,
        ReferencePosition pos,
        Optional<PileAllele> a1,
        Optional<PileAllele> a2) {
      super();
      this.pileup = pileup;
      Depth depth = pileup.getDepth();
      PileAllele ref = pos.getRefAllele();
      Optional<PileAllele> alt = pos.getAltAllele();

      this.id = id;
      rawDepth = depth.rawTotalDepth();
      refRawDepth = depth.allelicRawDepth(ref);
      altRawDepth = alt.transform(depth::allelicRawDepth);
      a1RawDepth = a1.transform(depth::allelicRawDepth).or(0);
      a2RawDepth = a2.transform(depth::allelicRawDepth).or(0);
      a_rawDepth = depth.allelicRawDepth(SNPAllele.A);
      t_rawDepth = depth.allelicRawDepth(SNPAllele.T);
      c_rawDepth = depth.allelicRawDepth(SNPAllele.C);
      g_rawDepth = depth.allelicRawDepth(SNPAllele.G);
      a1ClippedReads = a1.transform(pileup.getClippedReadCounts()::count).or(0);
      a2ClippedReads = a2.transform(pileup.getClippedReadCounts()::count).or(0);
      a1ApparentMismapReads = a1.transform(pileup.getApparentMismapReadCounts()::count).or(0);
      a2ApparentMismapReads = a2.transform(pileup.getApparentMismapReadCounts()::count).or(0);
      a1UnmappedMateReads = a1.transform(pileup.getUnmappedMateCounts()::count).or(0);
      a2UnmappedMateReads = a2.transform(pileup.getUnmappedMateCounts()::count).or(0);
      weightedDepth = depth.weightedTotalDepth();
      refWeightedDepth = depth.allelicWeightedDepth(ref);
      altWeightedDepth = alt.transform(depth::allelicWeightedDepth);
      a1WeightedDepth = a1.transform(depth::allelicWeightedDepth).or(0.0);
      a2WeightedDepth = a2.transform(depth::allelicWeightedDepth).or(0.0);
      a_weightedDepth = depth.allelicWeightedDepth(SNPAllele.A);
      t_weightedDepth = depth.allelicWeightedDepth(SNPAllele.T);
      c_weightedDepth = depth.allelicWeightedDepth(SNPAllele.C);
      g_weightedDepth = depth.allelicWeightedDepth(SNPAllele.G);
    }

    /** @return the id */
    public String getId() {
      return id;
    }

    /** @return the depth */
    public Depth getDepth() {
      return pileup.getDepth();
    }

    /** @return the pileup */
    public Pileup getPileup() {
      return pileup;
    }
  }

  private static final int MIN_PARENTAL_DEPTH = 10;

  private static final int MAX_PARENTAL_ALLELIC_DEPTH = 1;

  private static final double MAX_PARENTAL_ALLELIC_FRAC = 0.5;

  private static final String NO_NON_SUPERNOVO_REASON = ".";

  private static final String SNP_EFF_JAR = "/home/pankrat2/public/bin/snpEff/snpEff.jar";
  private static final String SNPEFF_GENOME = "hg19";
  private static final String SNPEFF_ANN_FIELD = "ANN";
  private static final String SNPEFF_ANN_DELIM = "[\\s]?\\|[\\s]?";
  private static final String SNPEFF_GENE = "Gene_Name";
  private static final String SNPEFF_ANNO = "Annotation";
  private static final String SNPEFF_IMPACT = "Annotation_Impact";

  public final String chr;
  public final int position;
  public final PileAllele refAllele;
  public final Optional<PileAllele> altAllele;
  public final Optional<PileAllele> allele1;
  public final Optional<PileAllele> allele2;
  public boolean biallelicHeterozygote;
  public boolean deNovo;
  public boolean superNovo;
  public String nonSuperNovoReason;
  public final double meanHaplotypeConcordance;
  public final int overlappingReadsHetCount;
  public static final double MIN_HAPLOTYPE_CONCORDANCE = 0.75;
  public final int overlappingReadsDiscordantHetCount;
  public final int overlappingReadsAdjacentDeNovoCounts;
  public final int overlappingReadsIndependentDeNovoCount;
  public final int overlapingReadsThirdAlleleCount;
  public final Sample child;
  public final Sample p1;
  public final Sample p2;
  public String snpeffGene;
  public String snpeffAnnotation;
  public String snpeffImpact;

  private final ReferencePosition pos;
  private final HaplotypeEvaluator.Result hapResults;
  private final ImmutableList<Sample> parents;

  public DeNovoResult(
      ReferencePosition pos,
      HaplotypeEvaluator.Result hapResults,
      Sample child,
      Sample p1,
      Sample p2) {
    this.pos = pos;
    this.hapResults = hapResults;
    this.child = child;
    this.p1 = p1;
    this.p2 = p2;
    this.parents = ImmutableList.of(p1, p2);

    position = pos.getPosition();
    chr = pos.getContig();
    refAllele = pos.getRefAllele();
    altAllele = pos.getAltAllele();
    allele1 = child.getDepth().getA1();
    allele2 = child.getDepth().getA2();

    if (hapResults.getConcordances().isEmpty()) meanHaplotypeConcordance = 1.0;
    else
      meanHaplotypeConcordance =
          hapResults
              .getConcordances()
              .stream()
              .mapToDouble(Double::valueOf)
              .summaryStatistics()
              .getAverage();
    overlappingReadsHetCount = hapResults.getConcordances().size();
    overlappingReadsDiscordantHetCount =
        (int)
            hapResults
                .getConcordances()
                .stream()
                .mapToDouble(Double::valueOf)
                .filter(d -> d < MIN_HAPLOTYPE_CONCORDANCE)
                .count();
    overlappingReadsAdjacentDeNovoCounts = hapResults.getAdjacentDeNovos();
    overlappingReadsIndependentDeNovoCount = hapResults.getOtherDeNovos();
    overlapingReadsThirdAlleleCount = hapResults.getOtherTriallelics();
    setFlags();
  }

  private void readObject(ObjectInputStream ois) throws ClassNotFoundException, IOException {
    ois.defaultReadObject();
    setFlags();
  }

  private void setFlags() {
    biallelicHeterozygote = TrioEvaluator.looksBiallelic(child.getPileup());
    deNovo = TrioEvaluator.looksDenovo(child.getPileup(), p1.getPileup(), p2.getPileup());
    if (!biallelicHeterozygote) nonSuperNovoReason = "Not biallelic heterozygote";
    else if (!deNovo) nonSuperNovoReason = "Not denovo";
    else if (Math.min(p1.weightedDepth, p2.weightedDepth) < MIN_PARENTAL_DEPTH)
      nonSuperNovoReason = "Parental weighted depth < " + MIN_PARENTAL_DEPTH;
    else if (hapResults.getOtherDeNovos() != 0) nonSuperNovoReason = "Other denovos in region";
    else if (meanHaplotypeConcordance < MIN_HAPLOTYPE_CONCORDANCE)
      nonSuperNovoReason = "Haplotype Concordance < " + MIN_HAPLOTYPE_CONCORDANCE;
    else if (hapResults.getOtherTriallelics() != 0) nonSuperNovoReason = "Triallelics in region";
    else nonSuperNovoReason = NO_NON_SUPERNOVO_REASON;
    superNovo = nonSuperNovoReason.equals(NO_NON_SUPERNOVO_REASON);
  }

  public VariantContext generateVariantContext() {
    ImmutableList.Builder<Allele> allelesBuilder = ImmutableList.builder();
    allelesBuilder.add(Allele.create(refAllele.toString(), true));
    allelesBuilder.addAll(altAllele.transform(a -> Allele.create(a.toString(), false)).asSet());
    List<Allele> alleles = allelesBuilder.build();
    return new VariantContextBuilder()
        .source("SuperNovo")
        .chr(chr)
        .start(position)
        .computeEndFromAlleles(alleles, position)
        .alleles(alleles)
        .make();
  }

  private void setAnnos(VariantContext annotatedVC, Map<String, Integer> fieldIndices) {
    if (chr.equals(annotatedVC.getContig()) && position == annotatedVC.getStart()) {
      String[] annFields =
          annotatedVC.getAttributeAsString(SNPEFF_ANN_FIELD, "").split(SNPEFF_ANN_DELIM);
      if (annFields.length < 4) {
        App.LOG.warn("Skipping annotation of variant " + annotatedVC.toString());
      } else {
        //        snpeffGene = annFields[fieldIndices.get(SNPEFF_GENE)];
        //        snpeffAnnotation = annFields[fieldIndices.get(SNPEFF_ANNO)];
        //        snpeffImpact = annFields[fieldIndices.get(SNPEFF_IMPACT)];
        snpeffGene = annFields[3];
        snpeffAnnotation = annFields[1];
        snpeffImpact = annFields[2];
      }
    } else {
      throw new IllegalArgumentException(
          "Annotated VC not reading in correct order of DeNovoResults");
    }
  }

  public static void retrieveAnnos(List<DeNovoResult> deNovoResults) {
    App.LOG.log(Level.INFO, "Running SnpEff to annotate variants");
    List<String> cmd =
        ImmutableList.<String>builder()
            .add("java")
            .add("-Xmx8G")
            .add("-jar")
            .add(SNP_EFF_JAR)
            .add("eff")
            .add(SNPEFF_GENOME)
            .add("-")
            .build();

    try {
      Process snpEffProc = new ProcessBuilder(cmd).redirectError(Redirect.INHERIT).start();
      try (VariantContextWriter writer =
          new VariantContextWriterBuilder()
              .clearOptions()
              .clearIndexCreator()
              .setOutputVCFStream(new BufferedOutputStream(snpEffProc.getOutputStream()))
              .build()) {
        Thread readInThread =
            new Thread(
                () -> {
                  try (VCFIterator vcfIter =
                      new VCFIteratorBuilder()
                          .open(new BufferedInputStream(snpEffProc.getInputStream()))) {
                    App.LOG.info(vcfIter.getHeader().toString());
                    String annDescrip =
                        vcfIter.getHeader().getInfoHeaderLine(SNPEFF_ANN_FIELD).getDescription();
                    String[] annHeader =
                        annDescrip
                            .substring(annDescrip.indexOf('\'') + 1, annDescrip.lastIndexOf('\''))
                            .split(SNPEFF_ANN_DELIM);
                    App.LOG.info(Arrays.toString(annHeader));
                    Map<String, Integer> fieldIndices =
                        IntStream.range(0, annHeader.length)
                            .boxed()
                            .collect(ImmutableMap.toImmutableMap(i -> annHeader[i], i -> i));
                    int annotated = 0;
                    App.LOG.log(Level.INFO, "Reading in SnpEff annotations");
                    for (DeNovoResult deNovoResult : deNovoResults) {

                      deNovoResult.setAnnos(vcfIter.next(), fieldIndices);
                      if (annotated++ % 10000 == 0)
                        App.LOG.log(
                            Level.INFO,
                            "Annotated "
                                + annotated
                                + " variants (of "
                                + deNovoResults.size()
                                + ")");
                    }
                  } catch (IOException e) {
                    App.LOG.error(e);
                  }
                });
        readInThread.start();
        VCFHeader vcfHeader = new VCFHeader();
        writer.writeHeader(vcfHeader);
        App.LOG.log(Level.INFO, "Wrote VCF header to SnpEff, generating variants");
        int generated = 0;
        for (DeNovoResult deNovoResult : deNovoResults) {
          writer.add(deNovoResult.generateVariantContext());
          if (generated++ % 10000 == 0)
            App.LOG.log(
                Level.INFO,
                "Generated " + generated + " variants (of " + deNovoResults.size() + ")");
        }
      }

    } catch (IOException e) {
      App.LOG.error(e);
    }
  }
}
