package org.pankratzlab.supernovo;

import java.io.File;
import java.io.IOException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import com.google.common.base.Optional;
import picocli.CommandLine;
import picocli.CommandLine.ArgGroup;
import picocli.CommandLine.Option;

public class App implements Runnable {

  public static final Logger LOG = LogManager.getLogger(App.class);

  private static App APP = null;

  @Option(
      names = {"--vcf", "-v"},
      paramLabel = "VCF",
      description = "VCF with variants to query",
      required = true)
  private File vcf;

  @Option(
      names = {"--childBam", "--bam"},
      paramLabel = "BAM",
      description = "BAM of child",
      required = true)
  private File childBam;

  @Option(
      names = {"--childID", "--cID"},
      paramLabel = "ID",
      description = "Sample ID of child",
      required = true)
  private String childID;

  @ArgGroup(exclusive = true, multiplicity = "1")
  private TrioOrSolo trioOrSolo = new TrioOrSolo();

  public static class TrioOrSolo {
    @ArgGroup(exclusive = false, multiplicity = "0..1")
    public Parents parents;

    @Option(
        names = {"--solo", "-1"},
        description =
            "Run analysis on a single sample, otherwise must include parental information for Trio analysis")
    public boolean solo;
  }

  static class Parents {
    @Option(
        names = {"--parent1Bam", "--p1Bam"},
        paramLabel = "BAM",
        description = "BAM of parent 1",
        required = true)
    public File p1Bam;

    @Option(
        names = {"--parent1ID", "--p1ID"},
        paramLabel = "ID",
        description = "Sample ID of parent 1",
        required = true)
    public String p1ID;

    @Option(
        names = {"--parent2Bam", "--p2Bam"},
        paramLabel = "BAM",
        description = "BAM of parent 2",
        required = true)
    public File p2Bam;

    @Option(
        names = {"--parent2ID", "--p2ID"},
        paramLabel = "ID",
        description = "Sample ID of parent 2",
        required = true)
    public String p2ID;
  }

  @Option(
      names = {"--annovarDir", "-a"},
      paramLabel = "DIR",
      description = "Directory where annovar is located",
      required = false)
  private String annovarDir;

  @Option(
      names = {"--snpEff", "-s"},
      paramLabel = "JAR",
      description = "Path to snpeff jar",
      required = false)
  private String snpEffJar;

  @Option(
      names = {"--genome", "--snpEffGenome"},
      paramLabel = "GENOME",
      description = "Genome build argument for snpeff/annovar",
      required = true)
  private String snpEffGenome;

  @Option(
      names = {"--output", "-o"},
      paramLabel = "FILE",
      description = "Output file for parsed de novo variants",
      required = true)
  private File output;

  @Option(
      names = {"--vcfMaxParentAD"},
      paramLabel = "AD",
      description =
          "Maximum AD (Allelic Depth) value from VCF for De Novo Allele in a parent to evaluate variant. "
              + "Variants with parental AD above this value will be assumed inherited.",
      defaultValue = "4")
  private int vcfMaxParentAD;

  @Option(
      names = {"--minDepth"},
      paramLabel = "DEPTH",
      description = "Minimum weighted depth to consider calling a variant",
      defaultValue = "10")
  private int minDepth;

  @Option(
      names = {"--minAllelicDepth"},
      paramLabel = "DEPTH",
      description = "Minimum allelic depth to consider calling a variant",
      defaultValue = "4")
  private int minAllelicDepth;

  @Option(
      names = {"--minAllelicFrac"},
      paramLabel = "FRACTION",
      description = "Minimum allelic fraction to consider calling a variant",
      defaultValue = "0.1")
  private double minAllelicFrac;

  @Option(
      names = {"--minParentalDepth"},
      paramLabel = "DEPTH",
      description = "Minimum parental weighted depth to consider a de novo variant supernovo",
      defaultValue = "10")
  private int minParentalDepth;

  @Option(
      names = {"--minOtherDNAllelicDepth"},
      paramLabel = "DEPTH",
      description =
          "Minimum allelic depth to count a local variant in the other DNs in region count",
      defaultValue = "1.5")
  private double minOtherDNAllelicDepth;

  @Option(
      names = {"--minOtherDNAllelicDepthIndependent"},
      paramLabel = "DEPTH",
      description =
          "Minimum allelic depth to count a local variant in the other DNs in region count, if it fails minAllelicFrac",
      defaultValue = "3.0")
  private double minOtherDNAllelicDepthIndependent;

  @Option(
      names = {"--maxMiscallFrac"},
      paramLabel = "RATIO",
      description =
          "Maximum allelic fraction in parents to consider as miscalled bases. "
              + "Variants with parental allelic fraction for the putatively de novo allele above this value will be assumed inherited.",
      defaultValue = "0.05")
  private double maxMiscallFrac;

  @Option(
      names = {"--maxMiscallWeight"},
      paramLabel = "WEIGHT",
      description =
          "Maximum weighted depth in parents to consider as miscalled bases. "
              + "Variants with parental weighted depth for the putatively de novo allele  above this value will be assumed inherited.",
      defaultValue = "1.0")
  private double maxMiscallWeight;

  @Option(
      names = {"--minHaplotypeConcordance"},
      paramLabel = "CONORDANCE",
      description =
          "Minimum concordance with inherited variant haplotypes in the region. "
              + "Concordance is defined as the ratio of calls for alleles of local variants that "
              + "are concordant with the most likely haplotype of the De Novo variant to the "
              + "total number of calls for local variants on overlapping reads",
      defaultValue = "0.75")
  private double minHaplotypeConcordance;

  @Option(
      names = {"--haplotypeSearchDistance"},
      paramLabel = "DIST",
      description =
          "Distance to search up and down stream for variants on reads that overlap a potential De Novo variant. "
              + "Typically, this should be set to your read length (or greater for a more exhaustive search.)",
      defaultValue = "150")
  private int haplotypeSearchDistance;

  // Force the currently running App to always be statically available
  private App() {
    APP = this;
  }

  public static App getInstance() {
    return APP;
  }

  public static void main(String[] args) {
    new CommandLine(new App()).execute(args);
  }

  @Override
  public void run() {
    try {
      final Evaluator eval;
      if (trioOrSolo.solo) {
        eval = new VariantEvaluator(childBam);
      } else {
        eval = new TrioEvaluator(childBam, trioOrSolo.parents.p1Bam, trioOrSolo.parents.p2Bam);
      }
      eval.run(vcf, getOutput());
    } catch (IOException | ClassNotFoundException e) {
      LOG.error("An IO error was encountered", e);
    }
  }

  /** @return the childBam */
  File getChildBam() {
    return childBam;
  }

  /** @return the childID */
  String getChildID() {
    return childID;
  }

  /** @return the p1Bam */
  File getP1Bam() {
    return trioOrSolo.parents.p1Bam;
  }

  /** @return the p1ID */
  String getP1ID() {
    return trioOrSolo.parents.p1ID;
  }

  /** @return the p2Bam */
  File getP2Bam() {
    return trioOrSolo.parents.p2Bam;
  }

  /** @return the p2ID */
  String getP2ID() {
    return trioOrSolo.parents.p2ID;
  }

  /** @return the output */
  File getOutput() {
    return output;
  }

  /** @return the vcfMaxParentAD */
  int getVcfMaxParentAD() {
    return vcfMaxParentAD;
  }

  /** @return the minDepth */
  int getMinDepth() {
    return minDepth;
  }

  /** @return the minAllelicDepth */
  int getMinAllelicDepth() {
    return minAllelicDepth;
  }

  /** @return the minAllelicFrac */
  double getMinAllelicFrac() {
    return minAllelicFrac;
  }

  /** @return the maxMiscallFrac */
  double getMaxMiscallFrac() {
    return maxMiscallFrac;
  }

  /** @return the maxMiscallWeight */
  double getMaxMiscallWeight() {
    return maxMiscallWeight;
  }

  /** @return the minParentalDepth */
  public int getMinParentalDepth() {
    return minParentalDepth;
  }

  /** @return the minOtherDNAllelicDepth */
  double getMinOtherDNAllelicDepth() {
    return minOtherDNAllelicDepth;
  }

  /** @return the minOtherDNAllelicDepthIndependent */
  double getMinOtherDNAllelicDepthIndependent() {
    return minOtherDNAllelicDepthIndependent;
  }

  /** @return the minHaplotypeConcordance */
  public double getMinHaplotypeConcordance() {
    return minHaplotypeConcordance;
  }

  /** @return the annovarDir */
  public Optional<String> getAnnovarDir() {
    return Optional.fromNullable(annovarDir);
  }

  /** @return the snpEffJar */
  public Optional<String> getSnpEffJar() {
    return Optional.fromNullable(snpEffJar);
  }

  /** @return the snpEffGenome */
  public String getGenome() {
    return snpEffGenome;
  }

  /** @return the haplotypeSearchDistance */
  public int getHaplotypeSearchDistance() {
    return haplotypeSearchDistance;
  }
}
