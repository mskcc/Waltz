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
package org.mskcc.juber.waltz.pileup;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.Frequency;
import org.mskcc.juber.genotype.Genotype;
import org.mskcc.juber.genotype.GenotypeEventType;
import org.mskcc.juber.genotype.GenotypeID;
import org.mskcc.juber.waltz.Waltz;
import org.mskcc.juber.waltz.pileup.processors.PileupProcessor;

import gnu.trove.list.array.TByteArrayList;
import gnu.trove.map.hash.TObjectIntHashMap;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;

/**
 * @author Juber Patel
 * 
 *         class that represents a pileup of reads in a very efficient way
 * 
 */
public class RegionPileup
{
	private IndexedFastaSequenceFile referenceFasta;
	private int readPairMismatchPolicy;
	private byte[] referenceBases;
	private Interval interval;
	// the last valid position in the current pileup
	private int lastValidPositionIndex;
	private PositionPileup[] positions;
	private PositionPileup[] positionsWithoutDuplicates;
	private int insertMin;
	private int insertMax;
	private boolean mateUnmapped;
	private boolean mateDistanceUnexpected;
	private byte[] readBases;
	private byte[] baseQualities;
	private int readIndex;
	private boolean duplicate = false;
	/**
	 * map of fragment name to fragment span
	 */
	private Map<String, Fragment> fragments;
	private int pileupIndex;
	/**
	 * the first valid position in the pileup for the current read
	 */
	private int validPileupStart;

	/**
	 * Genotypes and fragments that have those genotypes.
	 */
	private Map<GenotypeID, Set<String>> genotypes;

	public RegionPileup(IndexedFastaSequenceFile referenceFasta,
			int maxIntervalLength, int insertMin, int insertMax,
			int readPairMismatchPolicy)
	{
		this.referenceFasta = referenceFasta;
		this.insertMin = insertMin;
		this.insertMax = insertMax;
		this.readPairMismatchPolicy = readPairMismatchPolicy;
		this.positions = new PositionPileup[maxIntervalLength];
		this.positionsWithoutDuplicates = new PositionPileup[maxIntervalLength];

		// initialize position pileups
		for (int i = 0; i < positions.length; i++)
		{
			positions[i] = new PositionPileup();
			positionsWithoutDuplicates[i] = new PositionPileup();
		}

		genotypes = new HashMap<GenotypeID, Set<String>>();
		fragments = new HashMap<String, Fragment>();
	}

	/**
	 * must call this method before using the region pileup
	 * 
	 * @param interval
	 */
	public void prepFor(Interval interval)
	{
		int end = interval.getEnd();
		this.lastValidPositionIndex = end - interval.getStart();

		referenceBases = referenceFasta.getSubsequenceAt(interval.getContig(),
				interval.getStart(), end).getBases();
		this.interval = interval;

		// clean the pileup for reuse
		for (int i = 0; i <= lastValidPositionIndex; i++)
		{
			positions[i].reset(referenceBases[i]);
			positionsWithoutDuplicates[i].reset(referenceBases[i]);
		}

		genotypes.clear();
		fragments.clear();
	}

	public void addRecord(SAMRecord record)
	{
		// set whether this read is marked as a duplicate
		duplicate = record.getDuplicateReadFlag();

		if (record.getReadPairedFlag())
		{
			setMateInfo(record);
		}

		// the overlap with the mate that should be ignored in processing
		int[] mateOverlap = null;

		readBases = record.getReadBases();
		baseQualities = record.getBaseQualities();

		processRecord(record, mateOverlap, baseQualities);
	}

	/**
	 * process the alignment record
	 * 
	 * 
	 * @param record
	 * @param mateOverlap
	 * @param baseQualities
	 */
	private void processRecord(SAMRecord record, int[] mateOverlap,
			byte[] baseQualities)
	{
		// example of CIGAR string: 146M3I4M5S
		// Assuming that CIGAR will only have these operators: M I D S H

		// TODO HANDLE NEGATIVE STRAND!!!

		MappedRead mappedRead = addToFragment(record);

		// currently not using it by making validPileupStart = 0.
		// the part of current read that overlaps with this interval also
		// completely overlaps with the paired read and has already been
		// processed.
		if (validPileupStart > lastValidPositionIndex)
		{
			return;
		}

		// tracks the current position in the read
		// this has to start at 0 because we are parsing the CIGAR string and
		// going through it.
		readIndex = 0;
		int operatorLength = 0;
		String fragmentName = record.getReadName();
		Cigar cigar = record.getCigar();
		List<CigarElement> elements = cigar.getCigarElements();

		for (int i = 0; i < elements.size(); i++)
		{
			CigarElement e = elements.get(i);
			CigarOperator operator = e.getOperator();
			operatorLength = e.getLength();

			if (operator.equals(CigarOperator.MATCH_OR_MISMATCH))
			{
				MatchMismatchRecord matchMismatchRecord = null;

				// add the bases
				for (int j = 0; j < operatorLength; j++)
				{
					if (pileupIndex >= validPileupStart
							&& pileupIndex <= lastValidPositionIndex)
					{
						// add to the unrolled match-mismatch record
						if (matchMismatchRecord == null)
						{
							matchMismatchRecord = new MatchMismatchRecord(
									interval.getContig(),
									interval.getStart() + pileupIndex,
									operatorLength);
						}

						matchMismatchRecord.add(
								positions[pileupIndex].getRefBase(),
								readBases[readIndex]);

						positions[pileupIndex].addBase(
								(char) readBases[readIndex], fragmentName,
								readPairMismatchPolicy);
						if (!duplicate)
						{
							positionsWithoutDuplicates[pileupIndex].addBase(
									(char) readBases[readIndex], fragmentName,
									readPairMismatchPolicy);
						}
					}

					// increment both pileup index and read index
					pileupIndex++;
					readIndex++;
				}

				if (matchMismatchRecord != null)
				{
					matchMismatchRecord.recordSubstitutions(mappedRead);
				}
			}
			else if (operator.equals(CigarOperator.INSERTION))
			{
				// TODO replace this boundary check with proper tracking of
				// CIGAR and quitting when it goes out of the region
				if (pileupIndex > validPileupStart
						&& pileupIndex <= lastValidPositionIndex)
				{
					positions[pileupIndex].addInsertion(operatorLength);
					if (!duplicate)
					{
						positionsWithoutDuplicates[pileupIndex]
								.addInsertion(operatorLength);
					}

					// add insertion to special genotypes map
					// make genotype id
					int precedingGenomicPosition = interval.getStart()
							+ (pileupIndex - 1);
					byte[] ref = new byte[] { referenceBases[pileupIndex - 1] };
					byte[] alt = new byte[operatorLength + 1];
					alt[0] = ref[0];
					System.arraycopy(readBases, readIndex, alt, 1,
							operatorLength);
					// copy(readBases, readIndex - 1,
					// readIndex + operatorLength);
					GenotypeID genotypeID = new GenotypeID(
							GenotypeEventType.INSERTION, interval.getContig(),
							precedingGenomicPosition, ref, alt);
					// add
					// addGenotype(genotypeID, fragmentName);
					mappedRead.addGenotype(genotypeID);
				}

				// increment readIndex but not PileupIndex
				readIndex += operatorLength;
			}
			else if (operator.equals(CigarOperator.DELETION))
			{
				// add deletion to the genotypes
				if (pileupIndex > validPileupStart
						&& pileupIndex <= lastValidPositionIndex)
				{
					// make genotype id
					int precedingGenomicPosition = interval.getStart()
							+ (pileupIndex - 1);
					byte[] alt = new byte[] { referenceBases[pileupIndex - 1] };
					byte[] ref = referenceFasta
							.getSubsequenceAt(interval.getContig(),
									precedingGenomicPosition,
									precedingGenomicPosition + operatorLength)
							.getBases();
					// byte[] ref = Arrays.copyOfRange(referenceBases,
					// pileupIndex - 1, pileupIndex + operatorLength);

					GenotypeID genotypeID = new GenotypeID(
							GenotypeEventType.DELETION, interval.getContig(),
							precedingGenomicPosition, ref, alt);
					// add
					// addGenotype(genotypeID, fragmentName);
					mappedRead.addGenotype(genotypeID);
				}

				// add deletions to the pileup
				for (int j = 0; j < operatorLength; j++)
				{
					if (pileupIndex >= validPileupStart
							&& pileupIndex <= lastValidPositionIndex)
					{
						positions[pileupIndex].addBase('D', fragmentName,
								readPairMismatchPolicy);

						if (!duplicate)
						{
							positionsWithoutDuplicates[pileupIndex].addBase('D',
									fragmentName, readPairMismatchPolicy);
						}
					}

					// increment pileupIndex but don't increment readIndex
					pileupIndex++;
				}
			}
			else if (operator.equals(CigarOperator.SOFT_CLIP))
			{
				// soft clips do not count in alignment length, no change in
				// pileupIndex

				// soft clip at the beginning of alignment
				if (i == 0)
				{
					addClip(operatorLength, pileupIndex,
							baseQualities[operatorLength], true, false);
				}
				else
				{
					// soft clip at the end of the alignment
					addClip(operatorLength, pileupIndex - 1,
							baseQualities[baseQualities.length - operatorLength
									- 1],
							false, false);
				}

				// but increment the read index
				readIndex += e.getLength();
			}
			else if (operator.equals(CigarOperator.HARD_CLIP))
			{
				// hard clips do not count in alignment length, nor are the
				// clipped bases included in the read sequence
				// so no change in either pileupIndex or readIndex

				// hard clip at the beginning of alignment
				if (i == 0)
				{
					addClip(operatorLength, pileupIndex, baseQualities[0], true,
							true);
				}
				else
				{
					// hard clip at the end of the alignment
					addClip(operatorLength, pileupIndex - 1,
							baseQualities[baseQualities.length - 1], false,
							true);
				}
			}
			else
			{
				// increment both pileupIndex and readIndex properly
				pileupIndex += operatorLength;
				readIndex += operatorLength;
			}
		}

	}

	private void recordAndAdjustSpanDummy(SAMRecord record)
	{
		Fragment fragmentSpan = fragments.get(record.getReadName());
		if (fragmentSpan == null)
		{
			fragmentSpan = new Fragment(record.getDuplicateReadFlag());
			fragments.put(record.getReadName(), fragmentSpan);
		}

		MappedRead recordSpan = fragmentSpan.add(record);

		validPileupStart = 0;
		pileupIndex = record.getAlignmentStart() - interval.getStart();
	}

	private MappedRead addToFragment(SAMRecord record)
	{
		Fragment fragment = fragments.get(record.getReadName());
		if (fragment == null)
		{
			fragment = new Fragment(duplicate);
			fragments.put(record.getReadName(), fragment);
		}

		MappedRead mappedRead = fragment.add(record);

		// don't remove the already processed part, we need to process it again
		// and tally. Set proper pileupIndex
		validPileupStart = 0;
		pileupIndex = record.getAlignmentStart() - interval.getStart();

		return mappedRead;

		// // this code removes the already processed part properly
		// // current record is fully contained in previous records
		// // make sure this record is not processed
		// if (recordSpan == null)
		// {
		// validPileupStart = lastValidPositionIndex + 1;
		// return;
		// }
		//
		// // points to the current position in the pileup
		// pileupIndex = record.getAlignmentStart() - interval.getStart();
		// validPileupStart = recordSpan.start - interval.getStart();
		// // if the overlap is ending to the left of current pileup window, it
		// // doesn't affect us
		// if (validPileupStart < 0)
		// {
		// validPileupStart = 0;
		// }

	}

	/**
	 * 
	 * @param genotypeID
	 * 
	 *            Add fragment to appropriate genotype
	 */
	private void addGenotype(GenotypeID genotypeID, String fragmentName)
	{
		Set<String> fragments = genotypes.get(genotypeID);

		if (fragments == null)
		{
			fragments = new HashSet<String>();
			genotypes.put(genotypeID, fragments);
		}

		fragments.add(fragmentName);
	}

	private void addClip(int length, int pileupIndex, byte baseQuality,
			boolean before, boolean hard)
	{
		if (before)
		{
			// TODO replace this boundary check with proper tracking of
			// CIGAR and quitting when it goes out of the region

			if (pileupIndex >= 0 && pileupIndex <= lastValidPositionIndex)
			{
				positions[pileupIndex].addClipEnd(baseQuality, hard);
				if (!duplicate)
				{
					positionsWithoutDuplicates[pileupIndex]
							.addClipEnd(baseQuality, hard);
				}
			}

			// go to the start of the hard clip
			pileupIndex -= length;
		}
		else
		{
			if (pileupIndex >= 0 && pileupIndex <= lastValidPositionIndex)
			{
				positions[pileupIndex].addClipStart(baseQuality, hard);
				if (!duplicate)
				{
					positionsWithoutDuplicates[pileupIndex]
							.addClipStart(baseQuality, hard);
				}
			}

			// go to the actual start of the hard clip
			pileupIndex++;
		}

		// add clipping info to positions
		for (int j = 0; j < length && pileupIndex >= 0
				&& pileupIndex <= lastValidPositionIndex; j++)
		{
			positions[pileupIndex].addClip(hard);
			if (!duplicate)
			{
				positionsWithoutDuplicates[pileupIndex].addClip(hard);
			}

			pileupIndex++;
		}
	}

	private void setMateInfo(SAMRecord record)
	{
		mateUnmapped = record.getMateUnmappedFlag();

		if (mateUnmapped)
		{
			// no need to check for unexpected distance
			return;
		}
		// read not mapped in a proper pair orientation
		else if (!record.getProperPairFlag())
		{
			mateDistanceUnexpected = true;
		}
		else
		{
			int distance = Math.abs(record.getInferredInsertSize());

			if (distance < insertMin || distance > insertMax)
			{
				mateDistanceUnexpected = true;
			}
			else
			{
				mateDistanceUnexpected = false;
			}
		}
	}

	/**
	 * capture the characteristics of the read alignment that will be used
	 * in further processing/in building the pileup.
	 * 
	 * 
	 * collect read level info like mate unmapped, mate not mapped in right
	 * place, split read etc.
	 * 
	 * @param record
	 */
	private int[] collectReadLevelInfo(SAMRecord record)
	{
		readBases = record.getReadBases();
		baseQualities = record.getBaseQualities();

		// the overlap with the mate that should be ignored in processing
		int[] mateOverlap = mateOverlap(record);

		return mateOverlap;
	}

	private int[] mateOverlap(SAMRecord record)
	{
		// TODO Auto-generated method stub
		return null;
	}

	public void giveViewTo(PileupProcessor processor)
	{
		// finalize the base counts
		for (int i = 0; i <= lastValidPositionIndex; i++)
		{
			positions[i].computeCounts();
			positionsWithoutDuplicates[i].computeCounts();
		}

		// finalize genotype counts
		computeGenotypeCounts();

		RegionPileupView view = new RegionPileupView(referenceBases, interval,
				lastValidPositionIndex, positions, positionsWithoutDuplicates,
				genotypes, fragments, insertMin, insertMax);

		processor.setRegionPileupView(view);
	}

	/**
	 * compute genotype counts by bringing genotypes from each fragment and
	 * putting them together
	 * 
	 * @return
	 */
	private void computeGenotypeCounts()
	{
		genotypes.clear();

		for (String name : fragments.keySet())
		{
			Fragment fragment = fragments.get(name);

			Set<GenotypeID> fragmentGenotypes = fragment
					.getGenotypes(readPairMismatchPolicy);
			for (GenotypeID genotypeID : fragmentGenotypes)
			{
				Set<String> names = genotypes.get(genotypeID);
				if (names == null)
				{
					names = new HashSet<String>();
					genotypes.put(genotypeID, names);
				}

				names.add(name);
			}
		}

		
//		GenotypeID genotypeID = new GenotypeID(GenotypeEventType.INSERTION,
//				"11", 108203632, new byte[] { 71 }, new byte[] { 71, 84 });
//		Set<String> names = genotypes.get(genotypeID);
//		if(names == null)
//		{
//			return;
//		}
//		
//		for(String name : names)
//		{
//			System.out.println(name);
//		}
//		
//		int a = 5;
	}

	private class MatchMismatchRecord
	{
		public final String contig;
		public final int genomicStartPosition;
		private TByteArrayList refs;
		private TByteArrayList readBases;

		public MatchMismatchRecord(String contig, int genomicStartPosition,
				int initialLength)
		{
			this.contig = contig;
			this.genomicStartPosition = genomicStartPosition;
			refs = new TByteArrayList(initialLength);
			readBases = new TByteArrayList(initialLength);

		}

		public void add(byte ref, byte readBase)
		{
			refs.add(ref);
			readBases.add(readBase);
		}

		public void recordSNVs(MappedRead mappedRead)
		{
			int size = refs.size();
			for (int i = 0; i < size; i++)
			{
				byte ref = refs.get(i);
				byte readBase = readBases.get(i);

				// add substitution
				if (ref != readBase)
				{
					GenotypeID genotypeID = new GenotypeID(
							GenotypeEventType.SNV, contig,
							genomicStartPosition + i, new byte[] { ref },
							new byte[] { readBase });

					mappedRead.addGenotype(genotypeID);
				}
			}
		}

		public void recordSubstitutions(MappedRead mappedRead)
		{
			int size = refs.size();
			int i = 0;
			while (i < size)
			{
				// no processing needed for matches or N's
				if (readBases.get(i) == 'N' || refs.get(i) == readBases.get(i))
				{
					i++;
					continue;
				}

				// find the index of the next matching base
				int nextMatchIndex = size;
				int j = i + 1;
				for (; j < size; j++)
				{
					// break on match or N
					if (readBases.get(j) == 'N'
							|| refs.get(j) == readBases.get(j))
					{
						nextMatchIndex = j;
						break;
					}
				}

				// reached the end of record without finding a match
				if (j == size)
				{
					nextMatchIndex = j;
				}

				GenotypeID genotypeID = null;

				// SNV
				if (nextMatchIndex - i == 1)
				{
					genotypeID = new GenotypeID(GenotypeEventType.SNV, contig,
							genomicStartPosition + i,
							refs.toArray(i, nextMatchIndex - i),
							readBases.toArray(i, nextMatchIndex - i));
				}
				else
				{
					// MNV
					genotypeID = new GenotypeID(GenotypeEventType.MNV, contig,
							genomicStartPosition + i,
							refs.toArray(i, nextMatchIndex - i),
							readBases.toArray(i, nextMatchIndex - i));
				}

				mappedRead.addGenotype(genotypeID);

				i = nextMatchIndex;

			}

		}

	}
}
