package org.pankratzlab.supernovo.pileup;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.function.Consumer;
import java.util.stream.IntStream;
import org.pankratzlab.supernovo.GenomePosition;
import com.google.common.base.Functions;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;

public class PileupCache {

  private final File bam;
  private final LoadingCache<GenomePosition, Pileup> pileups;

  private static final CacheBuilder<Object, Object> PILEUP_CACHE_BUILDER =
      CacheBuilder.newBuilder().softValues();

  /** @param bam */
  public PileupCache(File bam) {
    super();
    this.bam = bam;
    this.pileups = PILEUP_CACHE_BUILDER.build(queryingPileupLoader(bam));
  }

  public Pileup get(GenomePosition pos) {
    return pileups.getUnchecked(pos);
  }

  public Map<GenomePosition, Pileup> getRange(GenomePosition start, GenomePosition stop) {
    return queryRange(start, stop);
  }

  private Map<GenomePosition, Pileup> queryRange(GenomePosition start, GenomePosition stop) {
    String contig = start.getContig();
    if (!stop.getContig().equals(contig))
      throw new IllegalArgumentException("Start and stop range must be on same chromosome");
    Map<GenomePosition, Pileup.Builder> pileupBuilders =
        IntStream.rangeClosed(start.getPosition(), stop.getPosition())
            .mapToObj(pos -> new GenomePosition(contig, pos))
            .collect(ImmutableMap.toImmutableMap(Functions.identity(), Pileup.Builder::new));
    Map<GenomePosition, Pileup> currentlyMapped = pileups.getAllPresent(pileupBuilders.keySet());
    if (currentlyMapped.keySet().containsAll(pileupBuilders.keySet())) return currentlyMapped;
    try (SamReader reader = SAMPositionQueryOverlap.SR_FACTORY.open(bam);
        CloseableIterator<SAMRecord> iterator =
            reader.queryOverlapping(start.getContig(), start.getPosition(), stop.getPosition())) {

      final Consumer<SAMRecord> addToAll =
          read -> pileupBuilders.values().forEach(b -> b.addRecord(read));
      iterator.forEachRemaining(addToAll);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    Map<GenomePosition, Pileup> queryPileups =
        ImmutableMap.copyOf(Maps.transformValues(pileupBuilders, Pileup.Builder::build));
    pileups.putAll(queryPileups);
    return queryPileups;
  }

  private static CacheLoader<GenomePosition, Pileup> queryingPileupLoader(final File bam) {
    return new CacheLoader<GenomePosition, Pileup>() {

      @Override
      public Pileup load(GenomePosition pos) {
        try (SAMPositionQueryOverlap spqo = new SAMPositionQueryOverlap(bam, pos)) {
          return new Pileup(spqo.getRecords(), pos);
        }
      }
    };
  }
}
