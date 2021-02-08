package org.pankratzlab.supernovo.output;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.lang.ProcessBuilder.Redirect;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import java.util.stream.Stream;
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
import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class DeNovoResult implements OutputFields, Serializable {

  /** */
  private static final long serialVersionUID = 1L;

  public static class Sample implements OutputFields, Serializable {

    /** */
    private static final long serialVersionUID = 2L;

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
    public final int a1LastReadPosition;
    public final int a2LastReadPosition;
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
      a1LastReadPosition = a1.transform(pileup.getLastPositionReadCounts()::count).or(0);
      a2LastReadPosition = a2.transform(pileup.getLastPositionReadCounts()::count).or(0);
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

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + a1ApparentMismapReads;
      result = prime * result + a1ClippedReads;
      result = prime * result + a1LastReadPosition;
      result = prime * result + a1RawDepth;
      result = prime * result + a1UnmappedMateReads;
      long temp;
      temp = Double.doubleToLongBits(a1WeightedDepth);
      result = prime * result + (int) (temp ^ (temp >>> 32));
      result = prime * result + a2ApparentMismapReads;
      result = prime * result + a2ClippedReads;
      result = prime * result + a2LastReadPosition;
      result = prime * result + a2RawDepth;
      result = prime * result + a2UnmappedMateReads;
      temp = Double.doubleToLongBits(a2WeightedDepth);
      result = prime * result + (int) (temp ^ (temp >>> 32));
      result = prime * result + a_rawDepth;
      temp = Double.doubleToLongBits(a_weightedDepth);
      result = prime * result + (int) (temp ^ (temp >>> 32));
      result = prime * result + ((altRawDepth == null) ? 0 : altRawDepth.hashCode());
      result = prime * result + ((altWeightedDepth == null) ? 0 : altWeightedDepth.hashCode());
      result = prime * result + c_rawDepth;
      temp = Double.doubleToLongBits(c_weightedDepth);
      result = prime * result + (int) (temp ^ (temp >>> 32));
      result = prime * result + g_rawDepth;
      temp = Double.doubleToLongBits(g_weightedDepth);
      result = prime * result + (int) (temp ^ (temp >>> 32));
      result = prime * result + ((id == null) ? 0 : id.hashCode());
      result = prime * result + ((pileup == null) ? 0 : pileup.hashCode());
      result = prime * result + rawDepth;
      result = prime * result + refRawDepth;
      temp = Double.doubleToLongBits(refWeightedDepth);
      result = prime * result + (int) (temp ^ (temp >>> 32));
      result = prime * result + t_rawDepth;
      temp = Double.doubleToLongBits(t_weightedDepth);
      result = prime * result + (int) (temp ^ (temp >>> 32));
      temp = Double.doubleToLongBits(weightedDepth);
      result = prime * result + (int) (temp ^ (temp >>> 32));
      return result;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (!(obj instanceof Sample)) return false;
      Sample other = (Sample) obj;
      if (a1ApparentMismapReads != other.a1ApparentMismapReads) return false;
      if (a1ClippedReads != other.a1ClippedReads) return false;
      if (a1LastReadPosition != other.a1LastReadPosition) return false;
      if (a1RawDepth != other.a1RawDepth) return false;
      if (a1UnmappedMateReads != other.a1UnmappedMateReads) return false;
      if (Double.doubleToLongBits(a1WeightedDepth)
          != Double.doubleToLongBits(other.a1WeightedDepth)) return false;
      if (a2ApparentMismapReads != other.a2ApparentMismapReads) return false;
      if (a2ClippedReads != other.a2ClippedReads) return false;
      if (a2LastReadPosition != other.a2LastReadPosition) return false;
      if (a2RawDepth != other.a2RawDepth) return false;
      if (a2UnmappedMateReads != other.a2UnmappedMateReads) return false;
      if (Double.doubleToLongBits(a2WeightedDepth)
          != Double.doubleToLongBits(other.a2WeightedDepth)) return false;
      if (a_rawDepth != other.a_rawDepth) return false;
      if (Double.doubleToLongBits(a_weightedDepth)
          != Double.doubleToLongBits(other.a_weightedDepth)) return false;
      if (altRawDepth == null) {
        if (other.altRawDepth != null) return false;
      } else if (!altRawDepth.equals(other.altRawDepth)) return false;
      if (altWeightedDepth == null) {
        if (other.altWeightedDepth != null) return false;
      } else if (!altWeightedDepth.equals(other.altWeightedDepth)) return false;
      if (c_rawDepth != other.c_rawDepth) return false;
      if (Double.doubleToLongBits(c_weightedDepth)
          != Double.doubleToLongBits(other.c_weightedDepth)) return false;
      if (g_rawDepth != other.g_rawDepth) return false;
      if (Double.doubleToLongBits(g_weightedDepth)
          != Double.doubleToLongBits(other.g_weightedDepth)) return false;
      if (id == null) {
        if (other.id != null) return false;
      } else if (!id.equals(other.id)) return false;
      if (pileup == null) {
        if (other.pileup != null) return false;
      } else if (!pileup.equals(other.pileup)) return false;
      if (rawDepth != other.rawDepth) return false;
      if (refRawDepth != other.refRawDepth) return false;
      if (Double.doubleToLongBits(refWeightedDepth)
          != Double.doubleToLongBits(other.refWeightedDepth)) return false;
      if (t_rawDepth != other.t_rawDepth) return false;
      if (Double.doubleToLongBits(t_weightedDepth)
          != Double.doubleToLongBits(other.t_weightedDepth)) return false;
      if (Double.doubleToLongBits(weightedDepth) != Double.doubleToLongBits(other.weightedDepth))
        return false;
      return true;
    }
  }

  private static class AnnotationLibrary implements Serializable {

    private static final long serialVersionUID = 1L;

    private class Annotation implements OutputFields, Serializable {
      private static final long serialVersionUID = 1L;
      private final ImmutableMap<String, String> annotationValues;

      public Annotation(VariantContext vc) {
        annotationValues =
            annotationNames
                .stream()
                .collect(
                    ImmutableMap.toImmutableMap(
                        anno -> anno, anno -> vc.getAttributeAsString(anno, ".")));
      }

      @Override
      public Stream<String> fieldHeaders() {
        return annotationNames.stream();
      }

      @Override
      public Stream<String> fieldValues() {
        return fieldHeaders().map(annotationValues::get);
      }
    }

    private final ImmutableList<String> annotationNames;

    /** @param annotationNames */
    public AnnotationLibrary(ImmutableList<String> annotationNames) {
      super();
      this.annotationNames = annotationNames;
    }
  }

  private static final int MIN_PARENTAL_DEPTH = 10;

  private static final int MAX_PARENTAL_ALLELIC_DEPTH = 1;

  private static final double MAX_PARENTAL_ALLELIC_FRAC = 0.5;

  private static final String NO_NON_SUPERNOVO_REASON = ".";

  private static final String SNP_EFF_JAR = "/home/pankrat2/public/bin/snpEff/snpEff.jar";
  private static final String SNPEFF_ANN_FIELD = "ANN";
  private static final String SNPEFF_ANN_DELIM = "[\\s]?\\|[\\s]?";
  private static final String SNPEFF_GENE = "Gene_Name";
  private static final String SNPEFF_ANNO = "Annotation";
  private static final String SNPEFF_IMPACT = "Annotation_Impact";

  private static final String ANNOVAR_DIR = "/home/pankrat2/public/bin/ANNOVAR/annovar/";
  private static final String ANNOVAR_TABLE_COMMAND = "table_annovar.pl";

  public final String chr;
  public final int position;
  public transient String snpeffGene;

  public transient String snpeffAnnotation;

  public transient String snpeffImpact;
  public transient String snpeffHGVSc;
  public transient String snpeffHGVSp;

  public transient AnnotationLibrary.Annotation annovar;

  public final PileAllele refAllele;
  public final Optional<PileAllele> altAllele;
  public final Optional<PileAllele> allele1;
  public final Optional<PileAllele> allele2;
  public final Optional<PileAllele> dnAllele;
  public final Optional<Boolean> dnIsRef;
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
  private final ReferencePosition pos;
  private final HaplotypeEvaluator.Result hapResults;
  private final TrioEvaluator trioEvaluator;
  private final ImmutableList<Sample> parents;

  public DeNovoResult(
      ReferencePosition pos,
      HaplotypeEvaluator.Result hapResults,
      TrioEvaluator trioEvaluator,
      Sample child,
      Sample p1,
      Sample p2) {
    this.pos = pos;
    this.hapResults = hapResults;
    this.trioEvaluator = trioEvaluator;
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
    dnAllele = trioEvaluator.dnAllele(child.getPileup(), p1.getPileup(), p2.getPileup());
    dnIsRef = dnAllele.transform(refAllele::equals);

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
    biallelicHeterozygote = trioEvaluator.looksBiallelic(child.getPileup());
    deNovo = trioEvaluator.looksDenovo(child.getPileup(), p1.getPileup(), p2.getPileup());
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
    ImmutableList.Builder<Allele> allelesBuilder = ImmutableList.builderWithExpectedSize(2);
    Allele ref = Allele.create(refAllele.toString(), true);
    allelesBuilder.add(ref);
    altAllele
        .transform(a -> Allele.create(a.toString(), false))
        .toJavaUtil()
        .ifPresent(allelesBuilder::add);
    List<Allele> alleles = allelesBuilder.build();
    List<Allele> homRefAlleles = ImmutableList.of(ref, ref);
    return new VariantContextBuilder()
        .source("SuperNovo")
        .chr(chr)
        .start(position)
        .computeEndFromAlleles(alleles, position)
        .alleles(alleles)
        .genotypes(
            new GenotypeBuilder(child.getId(), alleles).make(),
            new GenotypeBuilder(p1.getId(), homRefAlleles).make(),
            new GenotypeBuilder(p2.getId(), homRefAlleles).make())
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

        snpeffHGVSc = annFields.length > 8 ? annFields[8] : ".";
        snpeffHGVSp = annFields.length > 9 ? annFields[9] : ".";
      }
    } else {
      throw new IllegalArgumentException(
          "Annotated VC not reading in correct order of DeNovoResults");
    }
  }

  public static void retrieveAnnos(
      Collection<DeNovoResult> deNovoResults,
      File vcfOutputRoot,
      SAMSequenceDictionary dictionary,
      String snpEffGenome) {
    File intermediateVCFOutput =
        new File(vcfOutputRoot.getParentFile(), "TEMP_" + vcfOutputRoot.getName() + ".vcf.gz");
    App.LOG.log(Level.INFO, "Running SnpEff to annotate variants");
    List<String> cmd =
        ImmutableList.<String>builder()
            .add("java")
            .add("-Xmx8G")
            .add("-jar")
            .add(SNP_EFF_JAR)
            .add("eff")
            .add(snpEffGenome)
            .add("-")
            .build();

    try (VariantContextWriter snpEffOutWriter =
        new VariantContextWriterBuilder().setOutputFile(intermediateVCFOutput).build()) {
      Process snpEffProc = new ProcessBuilder(cmd).redirectError(Redirect.INHERIT).start();
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
                  snpEffOutWriter.writeHeader(vcfIter.getHeader());
                  for (DeNovoResult deNovoResult : deNovoResults) {
                    VariantContext annotatedVC = vcfIter.next();
                    snpEffOutWriter.add(annotatedVC);
                    deNovoResult.setAnnos(annotatedVC, fieldIndices);
                    if (annotated++ % 10000 == 0)
                      App.LOG.log(
                          Level.INFO,
                          "Annotated " + annotated + " variants (of " + deNovoResults.size() + ")");
                  }
                  if (vcfIter.hasNext()) App.LOG.warn("VCF Iteration has not completed");
                } catch (IOException e) {
                  App.LOG.error(e);
                }
              });
      readInThread.start();
      try (VariantContextWriter snpEffPipeWriter =
          new VariantContextWriterBuilder()
              .clearOptions()
              .clearIndexCreator()
              .setOutputVCFStream(new BufferedOutputStream(snpEffProc.getOutputStream()))
              .build()) {
        ImmutableList<String> samples =
            deNovoResults
                .stream()
                .findFirst()
                .map(dnr -> ImmutableList.of(dnr.child.getId(), dnr.p1.getId(), dnr.p2.getId()))
                .orElseGet(ImmutableList::of);
        VCFHeader vcfHeader =
            new VCFHeader(ImmutableSet.of(VCFStandardHeaderLines.getFormatLine("GT")), samples);
        vcfHeader.setSequenceDictionary(dictionary);
        snpEffPipeWriter.writeHeader(vcfHeader);
        App.LOG.log(Level.INFO, "Wrote VCF header to SnpEff, generating variants");
        int generated = 0;
        for (DeNovoResult deNovoResult : deNovoResults) {
          snpEffPipeWriter.add(deNovoResult.generateVariantContext());
          if (generated++ % 10000 == 0)
            App.LOG.log(
                Level.INFO,
                "Generated " + generated + " variants (of " + deNovoResults.size() + ")");
        }
      }
      readInThread.join();
    } catch (IOException e) {
      App.LOG.error(e);
    } catch (InterruptedException e) {
      App.LOG.error(e);
      Thread.currentThread().interrupt();
    }
    App.LOG.info("Running annovar");
    List<String> annovarCmd =
        ImmutableList.of(
            ANNOVAR_DIR + ANNOVAR_TABLE_COMMAND,
            intermediateVCFOutput.getPath(),
            ANNOVAR_DIR + "humandb/",
            "-buildver",
            snpEffGenome,
            "-out",
            vcfOutputRoot.getPath(),
            "-remove",
            "-protocol",
            "refGene,cytoBand,exac03,avsnp147,dbnsfp30a",
            "-operation",
            "gx,r,f,f,f",
            "-nastring",
            ".",
            "-vcfinput",
            "--thread",
            String.valueOf(ForkJoinPool.getCommonPoolParallelism()));
    try {
      Process annovarProc =
          new ProcessBuilder(annovarCmd)
              .redirectError(Redirect.INHERIT)
              .redirectOutput(Redirect.INHERIT)
              .start();
      annovarProc.waitFor();
    } catch (IOException e) {
      App.LOG.error("An error occurreed while running annovar", e);
    } catch (InterruptedException e) {
      App.LOG.error(e);
      Thread.currentThread().interrupt();
    }
    App.LOG.info("Finished running annovar, reading annotations");
    try (VCFFileReader annovarVcfReader =
            new VCFFileReader(new File(vcfOutputRoot.getPath() + ".hg19_multianno.vcf"), false);
        CloseableIterator<VariantContext> annovarIter = annovarVcfReader.iterator()) {
      AnnotationLibrary annovarLibrary =
          new AnnotationLibrary(
              annovarVcfReader
                  .getFileHeader()
                  .getInfoHeaderLines()
                  .stream()
                  .filter(info -> info.getDescription().contains("ANNOVAR"))
                  .map(VCFInfoHeaderLine::getID)
                  .collect(ImmutableList.toImmutableList()));
      for (DeNovoResult result : deNovoResults) {
        if (annovarIter.hasNext()) {
          VariantContext vc = annovarIter.next();
          if (result.position == vc.getStart()) {
            result.annovar = annovarLibrary.new Annotation(vc);
          } else {
            App.LOG.error(
                "Annovar VCF and results are not matched, \n\tnext result: "
                    + result.pos
                    + "\n\tnext annovar: "
                    + vc.getContig()
                    + ":"
                    + vc.getStart());
          }

        } else {
          App.LOG.error("Annovar VCF ended before annotating all de novo results");
        }
      }
    }
  }

  /** @return the pos */
  public ReferencePosition getPos() {
    return pos;
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((allele1 == null) ? 0 : allele1.hashCode());
    result = prime * result + ((allele2 == null) ? 0 : allele2.hashCode());
    result = prime * result + ((altAllele == null) ? 0 : altAllele.hashCode());
    result = prime * result + (biallelicHeterozygote ? 1231 : 1237);
    result = prime * result + ((child == null) ? 0 : child.hashCode());
    result = prime * result + ((chr == null) ? 0 : chr.hashCode());
    result = prime * result + (deNovo ? 1231 : 1237);
    result = prime * result + ((dnAllele == null) ? 0 : dnAllele.hashCode());
    result = prime * result + ((dnIsRef == null) ? 0 : dnIsRef.hashCode());
    result = prime * result + ((hapResults == null) ? 0 : hapResults.hashCode());
    long temp;
    temp = Double.doubleToLongBits(meanHaplotypeConcordance);
    result = prime * result + (int) (temp ^ (temp >>> 32));
    result = prime * result + ((nonSuperNovoReason == null) ? 0 : nonSuperNovoReason.hashCode());
    result = prime * result + overlapingReadsThirdAlleleCount;
    result = prime * result + overlappingReadsAdjacentDeNovoCounts;
    result = prime * result + overlappingReadsDiscordantHetCount;
    result = prime * result + overlappingReadsHetCount;
    result = prime * result + overlappingReadsIndependentDeNovoCount;
    result = prime * result + ((p1 == null) ? 0 : p1.hashCode());
    result = prime * result + ((p2 == null) ? 0 : p2.hashCode());
    result = prime * result + ((parents == null) ? 0 : parents.hashCode());
    result = prime * result + ((pos == null) ? 0 : pos.hashCode());
    result = prime * result + position;
    result = prime * result + ((refAllele == null) ? 0 : refAllele.hashCode());
    result = prime * result + ((snpeffAnnotation == null) ? 0 : snpeffAnnotation.hashCode());
    result = prime * result + ((snpeffGene == null) ? 0 : snpeffGene.hashCode());
    result = prime * result + ((snpeffHGVSc == null) ? 0 : snpeffHGVSc.hashCode());
    result = prime * result + ((snpeffHGVSp == null) ? 0 : snpeffHGVSp.hashCode());
    result = prime * result + ((snpeffImpact == null) ? 0 : snpeffImpact.hashCode());
    result = prime * result + (superNovo ? 1231 : 1237);
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj instanceof DeNovoResult)) return false;
    DeNovoResult other = (DeNovoResult) obj;
    if (allele1 == null) {
      if (other.allele1 != null) return false;
    } else if (!allele1.equals(other.allele1)) return false;
    if (allele2 == null) {
      if (other.allele2 != null) return false;
    } else if (!allele2.equals(other.allele2)) return false;
    if (altAllele == null) {
      if (other.altAllele != null) return false;
    } else if (!altAllele.equals(other.altAllele)) return false;
    if (biallelicHeterozygote != other.biallelicHeterozygote) return false;
    if (child == null) {
      if (other.child != null) return false;
    } else if (!child.equals(other.child)) return false;
    if (chr == null) {
      if (other.chr != null) return false;
    } else if (!chr.equals(other.chr)) return false;
    if (deNovo != other.deNovo) return false;
    if (dnAllele == null) {
      if (other.dnAllele != null) return false;
    } else if (!dnAllele.equals(other.dnAllele)) return false;
    if (dnIsRef == null) {
      if (other.dnIsRef != null) return false;
    } else if (!dnIsRef.equals(other.dnIsRef)) return false;
    if (hapResults == null) {
      if (other.hapResults != null) return false;
    } else if (!hapResults.equals(other.hapResults)) return false;
    if (Double.doubleToLongBits(meanHaplotypeConcordance)
        != Double.doubleToLongBits(other.meanHaplotypeConcordance)) return false;
    if (nonSuperNovoReason == null) {
      if (other.nonSuperNovoReason != null) return false;
    } else if (!nonSuperNovoReason.equals(other.nonSuperNovoReason)) return false;
    if (overlapingReadsThirdAlleleCount != other.overlapingReadsThirdAlleleCount) return false;
    if (overlappingReadsAdjacentDeNovoCounts != other.overlappingReadsAdjacentDeNovoCounts)
      return false;
    if (overlappingReadsDiscordantHetCount != other.overlappingReadsDiscordantHetCount)
      return false;
    if (overlappingReadsHetCount != other.overlappingReadsHetCount) return false;
    if (overlappingReadsIndependentDeNovoCount != other.overlappingReadsIndependentDeNovoCount)
      return false;
    if (p1 == null) {
      if (other.p1 != null) return false;
    } else if (!p1.equals(other.p1)) return false;
    if (p2 == null) {
      if (other.p2 != null) return false;
    } else if (!p2.equals(other.p2)) return false;
    if (parents == null) {
      if (other.parents != null) return false;
    } else if (!parents.equals(other.parents)) return false;
    if (pos == null) {
      if (other.pos != null) return false;
    } else if (!pos.equals(other.pos)) return false;
    if (position != other.position) return false;
    if (refAllele == null) {
      if (other.refAllele != null) return false;
    } else if (!refAllele.equals(other.refAllele)) return false;
    if (snpeffAnnotation == null) {
      if (other.snpeffAnnotation != null) return false;
    } else if (!snpeffAnnotation.equals(other.snpeffAnnotation)) return false;
    if (snpeffGene == null) {
      if (other.snpeffGene != null) return false;
    } else if (!snpeffGene.equals(other.snpeffGene)) return false;
    if (snpeffHGVSc == null) {
      if (other.snpeffHGVSc != null) return false;
    } else if (!snpeffHGVSc.equals(other.snpeffHGVSc)) return false;
    if (snpeffHGVSp == null) {
      if (other.snpeffHGVSp != null) return false;
    } else if (!snpeffHGVSp.equals(other.snpeffHGVSp)) return false;
    if (snpeffImpact == null) {
      if (other.snpeffImpact != null) return false;
    } else if (!snpeffImpact.equals(other.snpeffImpact)) return false;
    if (superNovo != other.superNovo) return false;
    return true;
  }
}
