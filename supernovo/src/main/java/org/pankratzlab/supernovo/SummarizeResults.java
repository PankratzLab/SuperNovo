package org.pankratzlab.supernovo;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.common.base.Optional;
import com.google.common.collect.LinkedHashMultiset;
import com.google.common.collect.Multiset;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;

public class SummarizeResults {

  public static void main(String[] args) throws ClassNotFoundException, IOException {
    summarizeAllResults(
        new File(
            "/home/spectorl/shared/Project_Spector_Project_014/gvcf_Oct2019/supernovo/SN_All_Summarized.txt"),
        new File("/home/spectorl/shared/Project_Spector_Project_014/gvcf_Oct2019/supernovo/"),
        new File("/home/spectorl/shared/Project_UMGC_Project_126/gvcf_Oct2019/supernovo"));
  }

  private static void summarizeAllResults(File output, File... dirs)
      throws ClassNotFoundException, IOException {
    Table<String, String, Integer> entrySampleCounts = TreeBasedTable.create();
    for (File file :
        Arrays.stream(dirs)
            .map(
                dir -> dir.listFiles((fileDir, name) -> name.endsWith(TrioEvaluator.SER_EXTENSION)))
            .flatMap(Arrays::stream)
            .toArray(File[]::new)) {
      Multiset<String> counts = LinkedHashMultiset.create();
      TrioEvaluator.deserializeResults(file)
          .values()
          .stream()
          .filter(Optional::isPresent)
          .map(Optional::get)
          .filter(d -> d.superNovo)
          .forEach(
              dnr -> {
                counts.add("supernovo");
                if (dnr.snpeffGene != null) counts.add(dnr.snpeffGene + "_AnyImpact");
                if (dnr.snpeffImpact != null) counts.add(dnr.snpeffImpact);
                if ("MODERATE".equals(dnr.snpeffImpact) || "HIGH".equals(dnr.snpeffImpact)) {
                  counts.add("supernovo_damaging");
                  counts.add(dnr.snpeffGene);
                  if (!dnr.dnIsRef.or(Boolean.FALSE)) counts.add("supernovo_damaging_nonref");
                }
              });
      counts.forEachEntry((entry, count) -> entrySampleCounts.put(file.getName(), entry, count));
    }
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(output)))) {
      writer.print("sample\t");
      writer.println(entrySampleCounts.columnKeySet().stream().collect(Collectors.joining("\t")));
      Set<String> cols = entrySampleCounts.columnKeySet();
      for (String sample : entrySampleCounts.rowKeySet()) {
        writer.print(sample + "\t");
        writer.println(
            cols.stream()
                .map(e -> entrySampleCounts.get(sample, e))
                .map(Optional::fromNullable)
                .map(Optional::get)
                .map(i -> i.toString())
                .collect(Collectors.joining("\t")));
      }
    }
  }
}
