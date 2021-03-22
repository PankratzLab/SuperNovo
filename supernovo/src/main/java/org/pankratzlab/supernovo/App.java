package org.pankratzlab.supernovo;

import java.io.File;
import java.io.IOException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import picocli.CommandLine;
import picocli.CommandLine.Option;

public class App implements Runnable {

  public static final Logger LOG = LogManager.getLogger(App.class);

  private static App APP = null;

  @Option(
      names = {"--vcf", "-v"},
      paramLabel = "VCF",
      description = "VCF with variants to query for de novo mutations",
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

  @Option(
      names = {"--parent1Bam", "--p1Bam"},
      paramLabel = "BAM",
      description = "BAM of parent 1",
      required = true)
  private File p1Bam;

  @Option(
      names = {"--parent1ID", "--p1ID"},
      paramLabel = "ID",
      description = "Sample ID of parent 1",
      required = true)
  private String p1ID;

  @Option(
      names = {"--parent2Bam", "--p2Bam"},
      paramLabel = "BAM",
      description = "BAM of parent 2",
      required = true)
  private File p2Bam;

  @Option(
      names = {"--parent2ID", "--p2ID"},
      paramLabel = "ID",
      description = "Sample ID of parent 2",
      required = true)
  private String p2ID;

  @Option(
      names = {"--genome", "--snpEffGenome"},
      paramLabel = "GENOME",
      description = "Genome build argument for snpeff",
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
      names = {"--maxMiscallFrac"},
      paramLabel = "RATIO",
      description =
          "Maximum allelic fraction in parents to consider as miscalled bases. "
              + "Variants with parental allelic fraction for the putatively de novo allele above this value will be assumed inherited.",
      defaultValue = "0.05")
  private double maxMiscallFrac = 0.05;

  @Option(
      names = {"--maxMiscallWeight"},
      paramLabel = "WEIGHT",
      description =
          "Maximum weighted depth in parents to consider as miscalled bases. "
              + "Variants with parental weighted depth for the putatively de novo allele  above this value will be assumed inherited.",
      defaultValue = "1.0")
  private double maxMiscallWeight = 1.0;

  // Force the currently running App to always be statically available
  private App() {
    APP = this;
  }

  public static App getInstance() {
    return APP;
  }

  public static void main(String[] args) {
    CommandLine.run(new App(), args);
  }

  @Override
  public void run() {
    try {
      new TrioEvaluator().reportDeNovos(vcf, getOutput());
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
    return p1Bam;
  }

  /** @return the p1ID */
  String getP1ID() {
    return p1ID;
  }

  /** @return the p2Bam */
  File getP2Bam() {
    return p2Bam;
  }

  /** @return the p2ID */
  String getP2ID() {
    return p2ID;
  }

  /** @return the snpEffGenome */
  String getSnpEffGenome() {
    return snpEffGenome;
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
}
