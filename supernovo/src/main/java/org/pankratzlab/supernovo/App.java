package org.pankratzlab.supernovo;

import java.io.File;
import java.io.IOException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import picocli.CommandLine;
import picocli.CommandLine.Option;

public class App implements Runnable {

  public static final Logger LOG = LogManager.getLogger(App.class);

  @Option(
    names = {"--vcf", "-v"},
    paramLabel = "VCF",
    description = "VCF with variants to query for de novo mutations",
    required = true
  )
  private File vcf;

  @Option(
    names = {"--childBam", "--bam"},
    paramLabel = "BAM",
    description = "BAM of child",
    required = true
  )
  private File childBam;

  @Option(
    names = {"--childID", "--cID"},
    paramLabel = "ID",
    description = "Sample ID of child",
    required = true
  )
  private String childID;

  @Option(
    names = {"--parent1Bam", "--p1Bam"},
    paramLabel = "BAM",
    description = "BAM of parent 1",
    required = true
  )
  private File p1Bam;

  @Option(
    names = {"--parent1ID", "--p1ID"},
    paramLabel = "ID",
    description = "Sample ID of parent 1",
    required = true
  )
  private String p1ID;

  @Option(
    names = {"--parent2Bam", "--p2Bam"},
    paramLabel = "BAM",
    description = "BAM of parent 2",
    required = true
  )
  private File p2Bam;

  @Option(
    names = {"--parent2ID", "--p2ID"},
    paramLabel = "ID",
    description = "Sample ID of parent 2",
    required = true
  )
  private String p2ID;

  @Option(
    names = {"--genome", "--snpEffGenome"},
    paramLabel = "GENOME",
    description = "Genome build argument for snpeff",
    required = true
  )
  private String snpEffGenome;

  @Option(
    names = {"--output", "-o"},
    paramLabel = "FILE",
    description = "Output file for parsed de novo variants",
    required = true
  )
  private File output;

  public static void main(String[] args) {
    CommandLine.run(new App(), args);
  }

  @Override
  public void run() {
    try {
      new TrioEvaluator(childBam, childID, p1Bam, p1ID, p2Bam, p2ID, snpEffGenome)
          .reportDeNovos(vcf, output);
    } catch (IOException | ClassNotFoundException e) {
      LOG.error("An IO error was encountered", e);
    }
  }
}
