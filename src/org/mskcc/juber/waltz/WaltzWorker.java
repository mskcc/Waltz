/*******************************************************************************
 *
 * @author Juber Patel
 *
 *         Copyright (c) 2017 Innovation Lab, CMO, MSKCC.
 *
 *         This software was developed at the Innovation Lab, Center for
 *         Molecular Oncology,
 *         Memorial Sloan Kettering Cancer Center, New York, New York.
 *
 *         Licensed under the Apache License, Version 2.0 (the "License");
 *         you may not use this file except in compliance with the License.
 *         You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 *         Unless required by applicable law or agreed to in writing, software
 *         distributed under the License is distributed on an "AS IS" BASIS,
 *         WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 *         implied.
 *         See the License for the specific language governing permissions and
 *         limitations under the License.
 *******************************************************************************/
/**
 * 
 */
package org.mskcc.juber.waltz;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.mskcc.juber.alignment.filters.AlignmentFilter;
import org.mskcc.juber.genotype.GenotypeID;
import org.mskcc.juber.waltz.pileup.RegionPileup;
import org.mskcc.juber.waltz.pileup.processors.GenotypingProcessor;
import org.mskcc.juber.waltz.pileup.processors.PileupMetricsProcessor;
import org.mskcc.juber.waltz.pileup.processors.PileupProcessor;
import org.mskcc.juber.waltz.pileup.processors.SignatureFindingProcessor;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

/**
 * @author Juber Patel
 * 
 */
public class WaltzWorker
{
	private PileupProcessor processor;
	private SamReader reader;
	private IntervalList intervalList;
	private int insertMin;
	private int insertMax;
	private int processedReads;
	private int validReads;
	private WaltzOutput output;
	private IndexedFastaSequenceFile referenceFasta;
	private AlignmentFilter filter;
	private int readPairMismatchPolicy;

	public WaltzWorker(String module, int minimumMappingQuality, String bamFile,
			String bamIndexFile, File referenceFastaFile,
			IntervalList intervalList, int readPairMismatchPolicy,
			String moduleArgument, int[] insertSize, WaltzOutput output)
			throws IOException
	{
		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamInputResource resource = SamInputResource.of(new File(bamFile))
				.index(new File(bamIndexFile));
		this.reader = factory.open(resource);
		this.intervalList = intervalList;
		this.insertMin = insertSize[0];
		this.insertMax = insertSize[1];
		this.output = output;
		this.referenceFasta = new IndexedFastaSequenceFile(referenceFastaFile);
		this.readPairMismatchPolicy = readPairMismatchPolicy;
		setFilter(minimumMappingQuality);
		setProcessor(module, moduleArgument);
	}

	private void setFilter(int minimumMappingQuality)
	{
		filter = new AlignmentFilter(minimumMappingQuality);
	}

	private void setProcessor(String module, String moduleArgument)
			throws IOException
	{
		if (module.equals("PileupMetrics"))
		{
			processor = new PileupMetricsProcessor();
			output.enableForMetrics();
		}
		else if (module.equals("Genotyping"))
		{
			String lociFilePath = moduleArgument;
			processor = new GenotypingProcessor(new File(lociFilePath),
					referenceFasta);
			Set<Interval> genotypeIntervals = ((GenotypingProcessor) processor)
					.getGenotypesAsIntervals();
			adjustIntervalList(genotypeIntervals);

			// output.enableForMetrics();
			output.enableForGenotypes(
					((GenotypingProcessor) processor).getMafHeader());

		}
		else if (module.equals("SignatureFinding"))
		{
			processor = new SignatureFindingProcessor(moduleArgument);
			output.enableForSignatureFinding();
		}
		else
		{
			System.err
					.println("Pileup Processor Type not recognized: " + module);
			System.err.println("Aborting.");
			System.exit(1);
		}
	}

	/**
	 * Adjust intervals list for genotyping. Remove intervals that do not
	 * overlap with any of the mutations. Expand intervals appropriately that
	 * partially overlap with a mutation.
	 * 
	 * @param genotypes
	 */
	private void adjustIntervalList(Set<Interval> genotypeIntervals)
	{
		// pad 5 bases on each side of interval. This is to make sure splice
		// sites are captured even when only exon boundaries are given.
		intervalList = intervalList.padded(5, 5);
		intervalList = intervalList.uniqued();

		// make new interval list that contains the parts of each interval that
		// overlaps with all the mutations overlapping with that interval.
		IntervalList newIntervalList = new IntervalList(
				intervalList.getHeader());

		for (Interval interval : intervalList.getIntervals())
		{
			IntervalList t = new IntervalList(intervalList.getHeader());

			for (Interval genotypeInterval : genotypeIntervals)
			{
				if (interval.intersects(genotypeInterval))
				{
					t.add(genotypeInterval);
				}
			}

			// no given mutation overlaps with current interval
			if (t.size() == 0)
			{
				continue;
			}

			// create one spanning interval that spans all the genotypes
			// overlapping with the current interval
			t = t.sorted();
			List<Interval> list = t.getIntervals();
			String contig = list.get(0).getContig();
			int start = list.get(0).getStart();
			int end = list.get(list.size() - 1).getEnd();
			Interval newInterval = new Interval(contig, start, end, false,
					interval.getName());
			newInterval = newInterval.pad(5, 5);
			newIntervalList.add(newInterval);
		}

		// update interval list that will be used for processing
		intervalList = newIntervalList;
	}

	/**
	 * check if the interval contains/overlaps the genotype
	 * 
	 * @param interval
	 * @param genotype
	 * @return
	 */
	private boolean overlap(Interval interval1, Interval interval2)
	{
		if (!interval1.getContig().equals(interval2.getContig()))
		{
			return false;
		}

		// sort by start
		if (interval2.getStart() < interval1.getStart())
		{
			Interval t = interval1;
			interval1 = interval2;
			interval2 = t;
		}

		if (interval1.getEnd() < interval2.getStart())
		{
			return false;
		}

		return true;
	}

	public Boolean process() throws IOException
	{
		long start = System.currentTimeMillis();

		// find out the maximum interval size
		int maxIntervalLength = -1;
		Iterator<Interval> it = intervalList.iterator();
		while (it.hasNext())
		{
			Interval interval = it.next();
			int length = interval.getEnd() - interval.getStart() + 1;
			if (length > maxIntervalLength)
			{
				maxIntervalLength = length;
			}
		}

		RegionPileup pileup = new RegionPileup(referenceFasta,
				maxIntervalLength, insertMin, insertMax, readPairMismatchPolicy);

		// for each interval
		for (Interval interval : intervalList)
		{
			System.out.println(interval);

			pileup.prepFor(interval);
			SAMRecordIterator iterator = reader.queryOverlapping(
					interval.getContig(), interval.getStart(),
					interval.getEnd());

			// for each record
			while (iterator.hasNext())
			{
				SAMRecord record = iterator.next();
				processedReads++;
				if (!filter.isGoodAlignment(record))
				{
					continue;
				}

				validReads++;

				try
				{
					pileup.addRecord(record);
				}
				catch (Exception e)
				{
					System.err.println("Problem processing record:");
					System.err.println(record.getSAMString());
					System.err.println("Region: " + interval);
					e.printStackTrace();
					continue;
				}
			}

			iterator.close();

			pileup.giveViewTo(processor);
			processor.processRegion(output);
		}

		long time = System.currentTimeMillis() - start;

		System.out.println("Processed " + processedReads + " reads total in "
				+ (time * 1.0) / 1000 + " seconds");
		System.out.println(validReads + " valid reads");

		// clean up
		referenceFasta.close();
		pileup = null;
		this.reader = null;
		this.intervalList = null;

		return true;

	}
}
