package org.pankratzlab.supernovo.contam;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.stream.Collectors;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import picocli.CommandLine;
import picocli.CommandLine.Option;

public class App implements Runnable {

  public static final Logger LOG = LogManager.getLogger(App.class);

  @Option(
    names = {"--vcf", "-v"},
    paramLabel = "VCF",
    description = "gVCF with variants to query for allelic ratios",
    required = true
  )
  private File vcf;

  @Option(
    names = {"--bins", "-b"},
    paramLabel = "n",
    description = "Number of bins to use for allele rations (default: ${DEFAULT-VALUE})"
  )
  private int bins = 100;

  @Option(
    names = {"--minAltDepth", "-a"},
    paramLabel = "n",
    description = "Minimum alternative allele depth to count (default: ${DEFAULT-VALUE})"
  )
  private int minAltDepth = 3;

  @Option(
    names = {"--minDepth", "-d"},
    paramLabel = "n",
    description = "Minimum total depth to count (default: ${DEFAULT-VALUE})"
  )
  private int minDepth = 3;

  @Option(
    names = {"--threads", "-t"},
    paramLabel = "n",
    description = "Number of threads to use (default: ${DEFAULT-VALUE})"
  )
  private int threads = 1;

  @Option(
    names = {"--output", "-o"},
    paramLabel = "FILE",
    description = "Output file for parsed allele ratio bins",
    required = true
  )
  private File output;

  public static void main(String[] args) {
    CommandLine.run(new App(), args);
  }

  @Override
  public void run() {
    try (VCFFileReader vcfReader = new VCFFileReader(vcf);
        CloseableIterator<VariantContext> vcfIter = vcfReader.iterator()) {
      if (vcfReader.getFileHeader().getSampleNamesInOrder().size() == 1) {
        AlleleRatio alleleRatio =
            new AlleleRatio(AlleleRatio.calculateBins(bins), minAltDepth, minDepth);
        vcfIter.stream().map(vc -> vc.getGenotype(0)).forEach(alleleRatio::addVariant);
        try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(output)))) {
          String header =
              alleleRatio
                  .getBins()
                  .stream()
                  .map(Object::toString)
                  .collect(Collectors.joining("\t"));
          writer.println(header);
          String data =
              alleleRatio
                  .getBins()
                  .stream()
                  .mapToInt(alleleRatio.getAltFracBinCounts()::count)
                  .mapToObj(Integer::toString)
                  .collect(Collectors.joining("\t"));
          writer.println(data);
        } catch (IOException e) {
          LOG.error("IO error encountered", e);
        }
      } else {
        LOG.error(
            vcf.getPath()
                + " is not a gVCF, contains "
                + vcfReader.getFileHeader().getSampleNamesInOrder().size()
                + " total samples");
      }
    }
  }
}
