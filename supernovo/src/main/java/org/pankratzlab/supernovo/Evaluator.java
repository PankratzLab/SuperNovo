package org.pankratzlab.supernovo;

import java.io.File;
import java.io.IOException;

public interface Evaluator {

  void run(File vcf, File output) throws IOException, ClassNotFoundException;
}