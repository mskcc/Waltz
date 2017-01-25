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


package org.mskcc.juber.waltz.bam;

import org.mskcc.juber.waltz.Waltz;

/**
 * @author Juber Patel
 * 
 *         Parse MD TAG from bam file
 * 
 */
public class MDTagParserCode
{
	private int value;

	public MDTagParserCode()
	{

	}

	private void processMDTag(String MDTag, int pileupIndex, int[] mateOverlap)
	{
		// example of MD tag: 5C12T1A99A30
		// assume that MD tag will only have numbers (matches), nucleotides
		// (mismatches) and ^ followed by nucleotides (deletions)

		char c;
		int length = MDTag.length();
		int lastNumber = 0;
		int nonNumericStart = -1;
		boolean inDigits = false;
		readIndex = 0;

		// parse the MD tag char by char
		for (int i = 0; i < length; i++)
		{
			c = MDTag.charAt(i);

			if (Character.isDigit(c))
			{
				// first digit of a digits section
				if (!inDigits)
				{
					inDigits = true;
					// get ready for the next number
					lastNumber = 0;

					// end of a series of bases or symbols etc. (non-numeric)
					// place to parse bases, deletions etc.
					if (nonNumericStart != -1)
					{
						pileupIndex = processMDNonNumber(MDTag, nonNumericStart,
								i - 1, pileupIndex);
					}
				}

				lastNumber = lastNumber * 10 + Character.getNumericValue(c);
			}
			else
			{
				// start of non-numeric section
				if (inDigits)
				{
					inDigits = false;

					// process the last parsed number
					pileupIndex = processMDNumber(lastNumber, pileupIndex);

					nonNumericStart = i;
				}
			}

		}

		// process the very last number or non-number in the MD tag
		if (inDigits)
		{
			pileupIndex = processMDNumber(lastNumber, pileupIndex);
		}
		else
		{
			pileupIndex = processMDNonNumber(MDTag, nonNumericStart, length - 1,
					pileupIndex);
		}
	}

	/**
	 * basically, add refs
	 * 
	 * @param lastNumber
	 * @param pileupIndex
	 * @return
	 */
	private int processMDNumber(int lastNumber, int pileupIndex)
	{
		// add refs
		for (int j = 0; j < lastNumber; j++)
		{
			if (pileupIndex >= 0 && pileupIndex < lastValidPositionIndex)
			{
				positions[pileupIndex].addRef(mateUnmapped,
						mateDistanceUnexpected);
				if (!duplicate)
				{
					positionsWithoutDuplicates[pileupIndex].addRef(mateUnmapped,
							mateDistanceUnexpected);
				}
			}

			pileupIndex++;
		}

		return pileupIndex;
	}

	private int processMDNonNumber(String MDTag, int nonNumericStart,
			int nonNumericEnd, int pileupIndex)
	{
		// a deletion
		if (MDTag.charAt(nonNumericStart) == '^')
		{
			for (int j = nonNumericStart + 1; j <= nonNumericEnd; j++)
			{
				if (pileupIndex >= 0 && pileupIndex < lastValidPositionIndex)
				{

					positions[pileupIndex].addDeletion(mateUnmapped,
							mateDistanceUnexpected);
					if (!duplicate)
					{
						positionsWithoutDuplicates[pileupIndex].addDeletion(
								mateUnmapped, mateDistanceUnexpected);
					}
				}

				pileupIndex++;
			}

			String deletion = MDTag.substring(nonNumericStart + 1,
					nonNumericEnd + 1);
			recordDeletionSequence(deletion);

		}
		else
		{
			// one or more mismatches
			for (int j = nonNumericStart; j <= nonNumericEnd; j++)
			{
				if (pileupIndex >= 0 && pileupIndex < lastValidPositionIndex)
				{

					positions[pileupIndex].addNonRef(mateUnmapped,
							mateDistanceUnexpected);
					if (!duplicate)
					{
						positionsWithoutDuplicates[pileupIndex].addNonRef(
								mateUnmapped, mateDistanceUnexpected);
					}
				}

				recordSubstitution(MDTag.charAt(j));
				pileupIndex++;
			}

		}

		return pileupIndex;
	}

	private void recordSubstitution(char base)
	{
		String ref = new String(new char[] { base });

		Integer count = Waltz.substituions.get(ref);

		if (count == null)
		{
			count = 1;
		}
		else
		{
			count += 1;
		}

		Waltz.substituions.put(ref, count);
	}

	private void recordDeletionSequence(String deletion)
	{
		// record only homopolymers
		// if (deletion.length() == 1)
		// return;
		// for (int i = 1; i < deletion.length(); i++)
		// {
		// if (deletion.charAt(i) != deletion.charAt(i - 1))
		// {
		// return;
		// }
		// }

		Integer count = Waltz.deletions.get(deletion);

		if (count == null)
		{
			count = 1;
		}
		else
		{
			count += 1;
		}

		Waltz.deletions.put(deletion, count);
	}


}
