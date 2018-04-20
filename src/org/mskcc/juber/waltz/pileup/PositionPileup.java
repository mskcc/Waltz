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

/**
 * @author Juber Patel
 * 
 *         Pileup at a specific genomic position (e.g. chr 2 position 3000000)
 * 
 */
public class PositionPileup
{

	/**
	 * number of reads supporting various events, in the order Ref, non-Ref, A,
	 * C, G, T, N, soft clip, hard clip, insertion,
	 * deletion, skipped region from genome.
	 * 
	 * Insertion read count refers to insertion immediately before this
	 * position.
	 */
	// private int[] readCount;

	private byte refBase;
	// read counts supporting various base calls
	// counts for A, C, G, T and N for this position
	private int[] baseCounts;
	private int[] baseQualities;
	private int softClips;
	private int softClipStarts;
	private int softClipStartQuality;
	private int softClipEnds;
	private int softClipEndQuality;
	private int hardClips;
	private int hardClipStarts;
	private int hardClipStartQuality;
	private int hardClipEnds;
	private int hardClipEndQuality;
	private int insertions;
	private int deletions;
	private int deletionsQuality;
	private int skipped;
	private int mateUnmapped;
	private int mateDistanceUnexpected;
	private int readSplitPoint;

	public PositionPileup()
	{
		// set the values that are not set by default
		// this value should never be seen anywhere
		refBase = '?';
		baseCounts = new int[5];
		baseQualities = new int[5];
	}

	public void reset(byte refBase)
	{
		this.refBase = refBase;
		for (int i = 0; i < baseCounts.length; i++)
		{
			baseCounts[i] = 0;
			baseQualities[i] = 0;

		}

		softClips = 0;
		softClipStarts = 0;
		softClipEnds = 0;
		softClipStartQuality = 0;
		softClipEndQuality = 0;
		hardClips = 0;
		hardClipStarts = 0;
		hardClipStartQuality = 0;
		hardClipEnds = 0;
		hardClipEndQuality = 0;

		insertions = 0;
		deletions = 0;
		skipped = 0;
		mateUnmapped = 0;
		mateDistanceUnexpected = 0;
		readSplitPoint = 0;
	}

	public void addBase(byte base)
	{
		if (base == 'A' || base == 'a')
		{
			baseCounts[0]++;
		}
		else if (base == 'C' || base == 'c')
		{
			baseCounts[1]++;
		}
		else if (base == 'G' || base == 'g')
		{
			baseCounts[2]++;
		}
		else if (base == 'T' || base == 't')
		{
			baseCounts[3]++;
		}
		else
		{
			baseCounts[4]++;
		}
	}

	public int getCount(byte base)
	{
		if (base == 'A' || base == 'a')
		{
			return baseCounts[0];
		}
		else if (base == 'C' || base == 'c')
		{
			return baseCounts[1];
		}
		else if (base == 'G' || base == 'g')
		{
			return baseCounts[2];
		}
		else if (base == 'T' || base == 't')
		{
			return baseCounts[3];
		}
		else if (base == 'I')
		{
			return insertions;
		}
		else if (base == 'D')
		{
			return deletions;
		}
		else
		{
			return -1;
		}
	}

	public int getCoverage()
	{
		// TODO double check if deletions should be part of this value
		return baseCounts[0] + baseCounts[1] + baseCounts[2] + baseCounts[3]
				+ deletions;
	}

	public void addDeletion(boolean mateUnmapped,
			boolean mateDistanceUnexpected)
	{
		deletions++;

		updateMateInfo(mateUnmapped, mateDistanceUnexpected);
	}

	public void addInsertion(int length)
	{
		insertions++;
	}

	/**
	 * increment clipping read number
	 * 
	 * @param hard
	 *            hard or soft
	 */
	public void addClip(boolean hard)
	{
		if (hard)
		{
			hardClips++;
		}
		else
		{
			softClips++;
		}
	}

	/**
	 * 
	 * @param baseQuality
	 * @param hard
	 *            hard or soft
	 */
	public void addClipEnd(byte baseQuality, boolean hard)
	{
		if (hard)
		{
			hardClipEnds++;
			hardClipEndQuality += baseQuality;
		}
		else
		{
			softClipEnds++;
			softClipEndQuality += baseQuality;
		}
	}

	/**
	 * 
	 * @param baseQuality
	 * @param hard
	 *            hard or soft
	 */
	public void addClipStart(byte baseQuality, boolean hard)
	{
		if (hard)
		{
			hardClipStarts++;
			hardClipStartQuality += baseQuality;
		}
		else
		{
			softClipStarts++;
			softClipStartQuality += baseQuality;
		}
	}

	private void updateMateInfo(boolean mateUnmapped,
			boolean mateDistanceUnexpected)
	{
		if (mateUnmapped)
		{
			this.mateUnmapped++;
		}
		else if (mateDistanceUnexpected)
		{
			this.mateDistanceUnexpected++;
		}
	}

	public String toString()
	{
		StringBuilder builder = new StringBuilder();

		builder.append((char) refBase).append('\t');
		builder.append(getCoverage() + baseCounts[4]).append('\t');
		builder.append(baseCounts[0]).append('\t');
		builder.append(baseCounts[1]).append('\t');
		builder.append(baseCounts[2]).append('\t');
		builder.append(baseCounts[3]).append('\t');
		builder.append(insertions).append('\t');
		builder.append(deletions).append('\t');
		builder.append(softClipStarts).append('\t');
		builder.append(softClipEnds).append('\t');
		builder.append(hardClipStarts).append('\t');
		builder.append(hardClipEnds);

		return builder.toString();
	}

	public byte getRefBase()
	{
		return refBase;
	}

	public int getHardClipStarts()
	{
		return hardClipStarts;
	}

	public int getHardClipEnds()
	{
		return hardClipEnds;
	}
}
