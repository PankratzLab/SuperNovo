package org.pankratzlab.supernovo.contam;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import com.google.common.primitives.Ints;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import picocli.CommandLine;
import picocli.CommandLine.Option;

public class App implements Runnable {

  public static final Logger LOG = LogManager.getLogger(App.class);

  @Option(
      names = {"--vcf", "-v"},
      paramLabel = "VCF",
      description = "gVCF with variants to query for allelic ratios",
      required = true)
  private File vcf;

  @Option(
      names = {"--bins", "-b"},
      paramLabel = "n",
      description = "Number of bins to use for allele rations (default: ${DEFAULT-VALUE})")
  private int bins = 100;

  @Option(
      names = {"--minAllelicDepth", "-a"},
      paramLabel = "n",
      description = "Minimum allelic depth to count (default: ${DEFAULT-VALUE})")
  private int minAllelicDepth = 3;

  @Option(
      names = {"--minDepth", "-d"},
      paramLabel = "n",
      description = "Minimum total depth to count (default: ${DEFAULT-VALUE})")
  private int minDepth = 20;

  @Option(
      names = {"--threads", "-t"},
      paramLabel = "n",
      description = "Number of threads to use (default: ${DEFAULT-VALUE})")
  private int threads = 4;

  @Option(
      names = {"--output", "-o"},
      paramLabel = "FILE",
      description = "Output file for parsed allele ratio bins",
      required = true)
  private File output;

  private AlleleRatio alleleRatio;

  public static void main(String[] args) {
    CommandLine.run(new App(), args);
  }

  private static boolean isAutosome(String contig) {
    final String contigNumber;
    if (contig.startsWith("chr")) contigNumber = contig.substring(3);
    else contigNumber = contig;
    return Ints.tryParse(contigNumber) != null;
  }

  private void countContig(String contig) {
    try (VCFFileReader vcfReader = new VCFFileReader(vcf);
        CloseableIterator<VariantContext> vcfIter = vcfReader.query(contig, 1, Integer.MAX_VALUE)) {
      vcfIter.stream().map(vc -> vc.getGenotype(0)).forEach(alleleRatio::addVariant);
    }
  }

  @Override
  public void run() {
    alleleRatio = new AlleleRatio(AlleleRatio.calculateBins(bins), minAllelicDepth, minDepth);
    final VCFHeader vcfHeader;
    try (VCFFileReader vcfReader = new VCFFileReader(vcf)) {
      vcfHeader = vcfReader.getFileHeader();
    }
    if (vcfHeader.getSampleNamesInOrder().size() == 1) {
      try {
        new ForkJoinPool(threads)
            .submit(
                () ->
                    vcfHeader
                        .getContigLines()
                        .stream()
                        .map(VCFContigHeaderLine::getID)
                        .filter(App::isAutosome)
                        .forEach(this::countContig))
            .get();
      } catch (InterruptedException e1) {
        LOG.error(e1);
        Thread.currentThread().interrupt();
      } catch (ExecutionException e1) {
        LOG.error(e1);
      }
      try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(output)))) {
        String header =
            alleleRatio.getBins().stream().map(Object::toString).collect(Collectors.joining("\t"));
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
              + vcfHeader.getSampleNamesInOrder().size()
              + " total samples");
    }
  }
}
