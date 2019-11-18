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
package org.mskcc.juber.waltz;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.io.FilenameUtils;

import com.google.common.base.Splitter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

/**
 * 
 */

/**
 * @author Juber Patel
 * 
 */
public class Waltz
{
	// TODO compute genome size??
	// right now, hard coding the genome size that is used only for
	// reporting
	// excludes X, Y and MT
	public static final long genomeSize = 2881033286L;
	public static final int basePadding = 0;

	/**
	 * @param args
	 * @throws IOException
	 * @throws Exception
	 */
	public static void main(String[] args) throws IOException
	{
		// TODO Should this program be made single threaded?? That would be more
		// suitable for cluster run and it will also resolve some design issues.

		final String module = args[0];
		final int minimumMappingQuality = Integer.parseInt(args[1]);
		final String bamFile = args[2];
		final File referenceFastaFile = new File(args[3]);
		final File intervalsBedFile = new File(args[4]);
		int readPairMismatchPolicy = 0;
		if (args.length >= 6)
		{
			readPairMismatchPolicy = Integer.parseInt(args[5]);
		}

		String moduleArgument = null;
		if (args.length >= 7)
		{
			moduleArgument = args[6];
		}

		// must not see args[] beyond this point

		final String bamIndexFile = FilenameUtils.removeExtension(bamFile)
				+ ".bai";
		final String sampleName = FilenameUtils
				.removeExtension(new File(bamFile).getName());

		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamInputResource resource = SamInputResource.of(new File(bamFile))
				.index(new File(bamIndexFile));
		SamReader reader = factory.open(resource);
		SAMFileHeader header = reader.getFileHeader();
		// TODO implement this method properly! It is broken!
		int[] dummyInsertSize = estimateInsertSize(reader, 99.0);
		reader.close();

		// IntervalList[] inputIntervalLists = makeIntervalLists(chunkSize,
		// numberOfThreads, header);
		IntervalList[] inputIntervalLists = makeIntervalLists(intervalsBedFile,
				1, header);
		WaltzOutput output = new WaltzOutput(sampleName);

		long start = System.currentTimeMillis();

		IntervalList intervalList = inputIntervalLists[0];
		WaltzWorker worker = new WaltzWorker(module, minimumMappingQuality,
				bamFile, bamIndexFile, referenceFastaFile, intervalList,
				readPairMismatchPolicy, moduleArgument, dummyInsertSize,
				output);

		// execute the worker
		worker.process();

		output.close();
		long time = System.currentTimeMillis() - start;
		System.out.println(
				"Program finished in " + (time * 1.0) / 1000 + " seconds\n");

		/**
		 * writeVariationFrequencies(insertions, sampleName +
		 * "-insertions.txt");
		 * writeVariationFrequencies(deletions, sampleName + "-deletions.txt");
		 * writeVariationFrequencies(substituions,
		 * sampleName + "-substitutions.txt");
		 **/
	}

	/**
	 * make given number of interval lists from the given bed file
	 * 
	 * @param targetBed
	 * @param numberOfLists
	 * @param header
	 * @return
	 * @throws IOException
	 */
	private static IntervalList[] makeIntervalLists(File bedFile,
			int numberOfLists, SAMFileHeader header) throws IOException
	{
		boolean addChr = hasChr(header);

		IntervalList[] intervalLists = new IntervalList[numberOfLists];
		for (int i = 0; i < intervalLists.length; i++)
		{
			intervalLists[i] = new IntervalList(header);
		}

		Splitter tabSplitter = Splitter.on('\t');

		BufferedReader reader = new BufferedReader(new FileReader(bedFile));
		String line = null;
		int lineNumber = 1;
		while ((line = reader.readLine()) != null)
		{
			List<String> words = tabSplitter.splitToList(line);
			String chr = null;
			if (addChr)
			{
				chr = "chr" + words.get(0);
			}
			else
			{
				chr = words.get(0);
			}

			// adding padding
			int start = Integer.parseInt(words.get(1)) - basePadding;
			int end = Integer.parseInt(words.get(2)) + basePadding;

			Interval interval = new Interval(chr, start, end, false,
					words.get(4));
			// add intervals to interval lists equitably
			IntervalList intervalList = intervalLists[lineNumber
					% numberOfLists];
			intervalList.add(interval);
		}

		reader.close();
		return intervalLists;
	}

	private static boolean hasChr(SAMFileHeader header)
	{
		String name = header.getSequence(0).getSequenceName();
		if (name.startsWith("chr"))
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	private static int[] estimateInsertSize(SamReader reader,
			double percentCovered) throws IOException
	{
		// System.out.println("Estimating min and max insert sizes that contain
		// "
		// + percentCovered + "% of the inserts from sampled reads");

		// System.out.println("Returning dummy range: 124 593");
		return new int[] { 124, 593 };

		// TODO the interval/reads to be sampled should be decided more
		// systematically
		/*
		 * SAMRecordIterator iterator = reader.queryOverlapping("chr1",
		 * 100000000,
		 * 110000000);
		 * 
		 * List<Integer> insertsTemp = new ArrayList<Integer>();
		 * 
		 * while (iterator.hasNext())
		 * {
		 * SAMRecord record = iterator.next();
		 * 
		 * if (record.getReadUnmappedFlag() || record.getMateUnmappedFlag())
		 * continue;
		 * 
		 * int insert = record.getInferredInsertSize();
		 * 
		 * if (insert <= 0)
		 * continue;
		 * 
		 * insertsTemp.add(insert);
		 * }
		 * 
		 * iterator.close();
		 * 
		 * double[] inserts = new double[insertsTemp.size()];
		 * 
		 * for (int i = 0; i < inserts.length; i++)
		 * {
		 * inserts[i] = insertsTemp.get(i);
		 * }
		 * 
		 * insertsTemp = null;
		 * 
		 * Arrays.sort(inserts);
		 * 
		 * Percentile percentile = new Percentile();
		 * percentile.setData(inserts);
		 * double p = (100.0 - percentCovered) / 2;
		 * 
		 * double left = percentile.evaluate(p);
		 * double right = percentile.evaluate(100 - p);
		 * 
		 * System.out.println((int) left + "\t" + (int) right);
		 * 
		 * return new int[] { (int) left, (int) right };
		 * 
		 */
	}

	private static IntervalList[] makeIntervalLists(int chunkSize,
			int numberOfLists, SAMFileHeader header)
	{
		String[] sequencesOfInterest = getSeqsOfInterest();

		List<Interval> intervals = new ArrayList<Interval>();

		long totalBases = 0;

		// make intervals of specified chunk size
		for (int i = 0; i < sequencesOfInterest.length; i++)
		{
			String sequenceName = sequencesOfInterest[i];
			int size = header.getSequence(sequenceName).getSequenceLength();
			totalBases += size;

			int start = 1;
			int end = chunkSize;

			while (start <= size)
			{
				intervals.add(new Interval(sequenceName, start, end));
				start = end + 1;
				end = start + chunkSize - 1;
			}
		}

		// shuffle the intervals
		Collections.shuffle(intervals);

		// create interval lists
		IntervalList[] intervalLists = new IntervalList[numberOfLists];
		for (int i = 0; i < intervalLists.length; i++)
		{
			intervalLists[i] = new IntervalList(header);
		}

		// add intervals to interval lists equitably
		for (int i = 0; i < intervals.size(); i++)
		{
			Interval interval = intervals.get(i);

			intervalLists[i % numberOfLists].add(interval);
		}

		return intervalLists;
	}

	private static String[] getSeqsOfInterest()
	{
		// return new String[] { "gi|215104|gb|J02459.1|LAMCG" };

		String[] seqs = new String[22];

		for (int i = 0; i < 22; i++)
		{
			seqs[i] = "chr" + Integer.toString(i + 1);
		}

		return seqs;
	}

}
