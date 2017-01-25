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
package org.mskcc.juber.waltz.pileup;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.mskcc.juber.genotype.Genotype;
import org.mskcc.juber.genotype.GenotypeEventType;
import org.mskcc.juber.genotype.GenotypeID;
import org.mskcc.juber.waltz.Waltz;
import org.mskcc.juber.waltz.pileup.processors.PileupProcessor;

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
	 * holds special genotypes: multi-base events and insertions
	 * multi-base substitution not handled yet.
	 */
	private Map<GenotypeID, Genotype> specialGenotypes;

	public RegionPileup(IndexedFastaSequenceFile referenceFasta, int insertMin,
			int insertMax)
	{
		this.referenceFasta = referenceFasta;
		this.insertMin = insertMin;
		this.insertMax = insertMax;
		this.positions = new PositionPileup[Waltz.getMaxRegionLength()];
		this.positionsWithoutDuplicates = new PositionPileup[Waltz
				.getMaxRegionLength()];

		// initialize position pileups
		for (int i = 0; i < positions.length; i++)
		{
			positions[i] = new PositionPileup();
			positionsWithoutDuplicates[i] = new PositionPileup();
		}

		specialGenotypes = new HashMap<GenotypeID, Genotype>();
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

		specialGenotypes.clear();
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

		// processMDTag(MDTag, pileupIndex, mateOverlap);
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

		// tracks the current position in the read
		readIndex = 0;
		// points to the current position in the pileup
		int pileupIndex = record.getStart() - interval.getStart();
		int operatorLength = 0;
		Cigar cigar = record.getCigar();
		List<CigarElement> elements = cigar.getCigarElements();

		for (int i = 0; i < elements.size(); i++)
		{
			CigarElement e = elements.get(i);
			CigarOperator operator = e.getOperator();
			operatorLength = e.getLength();

			if (operator.equals(CigarOperator.MATCH_OR_MISMATCH))
			{
				// add the bases
				for (int j = 0; j < operatorLength; j++)
				{
					if (pileupIndex >= 0
							&& pileupIndex <= lastValidPositionIndex)
					{
						positions[pileupIndex].addBase(readBases[readIndex]);
						if (!duplicate)
						{
							positionsWithoutDuplicates[pileupIndex]
									.addBase(readBases[readIndex]);
						}
					}

					// increment both pileup index and read index
					pileupIndex++;
					readIndex++;
				}
			}
			else if (operator.equals(CigarOperator.INSERTION))
			{
				// TODO replace this boundary check with proper tracking of
				// CIGAR and quitting when it goes out of the region
				if (pileupIndex >= 1
						&& pileupIndex <= lastValidPositionIndex + 1)
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
					addSpecialGenotype(genotypeID);
				}

				// increment readIndex but not PileupIndex
				readIndex += operatorLength;
			}
			else if (operator.equals(CigarOperator.DELETION))
			{
				// add deletion to the special genotypes map iff it is a
				// multi-base deletion
				if (operatorLength > 1 && pileupIndex >= 1
						&& pileupIndex <= lastValidPositionIndex + 1)
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
					addSpecialGenotype(genotypeID);
				}

				// add deletions to the pileup
				for (int j = 0; j < operatorLength; j++)
				{
					if (pileupIndex >= 0
							&& pileupIndex <= lastValidPositionIndex)
					{
						positions[pileupIndex].addDeletion(mateUnmapped,
								mateDistanceUnexpected);
						if (!duplicate)
						{
							positionsWithoutDuplicates[pileupIndex].addDeletion(
									mateUnmapped, mateDistanceUnexpected);
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

	/**
	 * 
	 * @param genotypeID
	 * 
	 *            increment the coverage count of the given genotype by 1
	 */
	private void addSpecialGenotype(GenotypeID genotypeID)
	{
		Genotype genotype = specialGenotypes.get(genotypeID);

		if (genotype == null)
		{
			genotype = new Genotype("FromPileup", genotypeID);
			specialGenotypes.put(genotypeID, genotype);
		}

		genotype.totalSupportingCoverage++;

		if (!duplicate)
		{
			genotype.uniqueSupportingCoverage++;
		}
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
		RegionPileupView view = new RegionPileupView(referenceBases, interval,
				lastValidPositionIndex, positions, positionsWithoutDuplicates,
				specialGenotypes, insertMin, insertMax);
		processor.setRegionPileupView(view);
	}
}
