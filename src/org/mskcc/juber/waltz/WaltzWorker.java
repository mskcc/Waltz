/*******************************************************************************
 *
 * @author Juber Patel
 *
 * Copyright (c) 2017 Innovation Lab, CMO, MSKCC.
 *
 * This software was developed at the Innovation Lab, Center for Molecular Oncology, 
 * Memorial Sloan Kettering Cancer Center, New York, New York.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *******************************************************************************/
/**
 * 
 */
package org.mskcc.juber.waltz;

import java.io.File;
import java.io.IOException;

import org.mskcc.juber.alignment.filters.AlignmentFilter;
import org.mskcc.juber.alignment.filters.BasicFilter;
import org.mskcc.juber.alignment.filters.MapQ20Filter;
import org.mskcc.juber.waltz.pileup.RegionPileup;
import org.mskcc.juber.waltz.pileup.processors.GenotypingProcessor;
import org.mskcc.juber.waltz.pileup.processors.MetricsProcessor;
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

	public WaltzWorker(String module, String filterType, String bamFile,
			String bamIndexFile, File referenceFastaFile,
			IntervalList intervalList, String moduleArgument, int[] insertSize,
			WaltzOutput output) throws IOException
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
		setFilter(filterType);
		setProcessor(module, moduleArgument);
	}

	private void setFilter(String filterType)
	{
		if (filterType.equals("BasicFilter"))
		{
			filter = new BasicFilter();
		}
		else if (filterType.equals("QualityFilter"))
		{
			filter = new MapQ20Filter();
		}
		else
		{
			System.err.println(
					"Alignment Filter Type not recognized: " + filterType);
			System.err.println("Aborting.");
			System.exit(1);
		}
	}

	private void setProcessor(String module, String moduleArgument)
			throws IOException
	{
		if (module.equals("Metrics"))
		{
			processor = new MetricsProcessor();
			output.enableForMetrics();
		}
		else if (module.equals("Genotyping"))
		{
			String lociFilePath = moduleArgument;
			processor = new GenotypingProcessor(new File(lociFilePath));
			output.enableForMetrics();
			output.enableForGenotypes();
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

	public Boolean process() throws IOException
	{
		long start = System.currentTimeMillis();
		RegionPileup pileup = new RegionPileup(referenceFasta, insertMin,
				insertMax);

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
					// System.exit(1);

				}
			}

			iterator.close();

			pileup.giveViewTo(processor);
			processor.processRegion(output);
		}

		// SamLocusIterator locusIterator = new SamLocusIterator(reader,
		// intervalList, true);
		//
		// // for each locus
		// while (locusIterator.hasNext())
		// {
		// SAMRecord record = null;
		// LocusInfo locusInfo = locusIterator.next();
		//
		// int positiveTotal = 0;
		// int negativeTotal = 0;
		// int[] positiveStrandBaseCounts = new int[4];
		// int[] negativeStrandBaseCounts = new int[4];
		// HashSet<String> seen = new HashSet<String>();
		//
		// // for each read covering the locus
		// for (RecordAndOffset read : locusInfo.getRecordAndPositions())
		// {
		// // ignore Ns
		// if (read.getReadBase() == 'N')
		// continue;
		//
		// // ignore bases below the quality threshold
		// if (read.getBaseQuality() < minBaseQuality)
		// continue;
		//
		// // choose the right strandwise accumulator
		// int[] rightArray = null;
		// record = read.getRecord();
		//
		// // count a read pair only once
		// if (seen.contains(record.getReadName()))
		// {
		// continue;
		// }
		//
		// seen.add(record.getReadName());
		//
		// if (record.getReadNegativeStrandFlag())
		// {
		// rightArray = negativeStrandBaseCounts;
		// negativeTotal++;
		// }
		// else
		// {
		// rightArray = positiveStrandBaseCounts;
		// positiveTotal++;
		// }
		//
		// rightArray[Utils.baseToInt(read.getReadBase()) - 1]++;
		//
		// }
		//
		// // no reads for current locus
		// if (record == null)
		// {
		// continue;
		// }

		// // add the pileup to the map
		// // Locus locus = new Locus(Utils.getChrNumber(record
		// // .getReferenceName()), locusInfo.getPosition());
		// // Pileup pileup = new Pileup(locusInfo.getSequenceName(),
		// // locusInfo.getPosition(), positiveStrandBaseCounts,
		// // negativeStrandBaseCounts, positiveTotal, negativeTotal);
		// // pileups.put(locus, pileup);
		// }
		//
		// locusIterator.close();
		// reader.close();

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
