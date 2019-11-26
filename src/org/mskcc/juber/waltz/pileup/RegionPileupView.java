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

import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.mskcc.juber.genotype.GenotypeID;

import htsjdk.samtools.util.Interval;

/**
 * @author Juber Patel
 *
 */
public class RegionPileupView
{
	public final byte[] referenceBases;
	public final Interval interval;
	// the last valid position in the current pileup
	public final int lastValidPositionIndex;
	public final PositionPileup[] positions;
	public final PositionPileup[] positionsWithoutDuplicates;
	/**
	 * holds special genotypes: multi-base events and insertions
	 * multi-base substitution not handled yet.
	 */
	public final Map<GenotypeID, Set<String>> genotypes;
	public final Map<String, FragmentSpan> fragmentSpans;
	public final int insertMin;
	public final int insertMax;

	public RegionPileupView(byte[] referenceBases, Interval interval,
			int lastValidPositionIndex, PositionPileup[] positions,
			PositionPileup[] positionsWithoutDuplicates,
			Map<GenotypeID, Set<String>> genotypes,
			Map<String, FragmentSpan> fragmentSpans, int insertMin,
			int insertMax)
	{
		this.referenceBases = referenceBases;
		this.interval = interval;
		this.lastValidPositionIndex = lastValidPositionIndex;
		this.positions = positions;
		this.positionsWithoutDuplicates = positionsWithoutDuplicates;
		this.genotypes = genotypes;
		this.fragmentSpans = fragmentSpans;
		this.insertMin = insertMin;
		this.insertMax = insertMax;
	}

	/**
	 * get the index of the given position in this pileup view
	 * 
	 * @param chr
	 * @param position
	 * @return index or -1 if it does not exist in this pileup
	 */
	public int getIndex(String chr, int position)
	{
		if (!interval.getContig().equals(chr))
		{
			return -1;
		}
		else if (interval.getStart() <= position
				&& interval.getEnd() >= position)
		{
			// calculate index
			return position - interval.getStart();
		}
		else
		{
			return -1;
		}
	}

	public boolean contains(GenotypeID genotypeID)
	{
		int pileupIndex = getIndex(genotypeID.contig, genotypeID.position);

		if (pileupIndex == -1)
		{
			return false;
		}

		pileupIndex = getIndex(genotypeID.contig, genotypeID.endPosition);

		if (pileupIndex == -1)
		{
			return false;
		}

		return true;
	}

	/**
	 * go through the position pileups and find the fragments that have valid
	 * coverage (i.e. no N's) over the entire given span
	 * 
	 * 
	 * @param span
	 * @return
	 */
	public Set<String> getValidSpanningFragments(GenotypeID genotypeID)
	{
		Set<String> spanningFragments = new HashSet<String>();

		int startIndex = getIndex(genotypeID.contig, genotypeID.position);
		int endIndex = getIndex(genotypeID.contig, genotypeID.endPosition);

		if (startIndex == -1 || endIndex == -1)
		{
			return spanningFragments;
		}

		// make the fragment set for the first position
		Map<String, Character> bases = positions[startIndex].getFragmentBases();
		for (String fragment : bases.keySet())
		{
			Character base = bases.get(fragment);

			// valid base
			if (base != null && base != PositionPileup.nChar)
			{
				spanningFragments.add(fragment);
			}
		}

		// remove those not present or not valid for other positions
		for (int i = startIndex + 1; i <= endIndex; i++)
		{
			// nothing to remove, short circuit
			if (spanningFragments.isEmpty())
			{
				return spanningFragments;
			}

			bases = positions[i].getFragmentBases();
			Iterator<String> iterator = spanningFragments.iterator();
			while (iterator.hasNext())
			{
				String fragment = iterator.next();
				Character base = bases.get(fragment);
				if (base == null || base == PositionPileup.nChar)
				{
					iterator.remove();
				}
			}
		}

		return spanningFragments;

	}

	public Set<String> getValidSpanningFragments(Set<GenotypeID> genotypeIDs)
	{
		Set<String> spanningFragments = null;
		for (GenotypeID genotypeID : genotypeIDs)
		{
			// first genotypeID
			if (spanningFragments == null)
			{
				spanningFragments = getValidSpanningFragments(genotypeID);
			}
			else if (spanningFragments.isEmpty())
			{
				// nothing to remove, return
				return spanningFragments;
			}
			else
			{
				Set<String> next = getValidSpanningFragments(genotypeID);
				// do intersection
				spanningFragments.retainAll(next);
			}
		}

		if (spanningFragments == null)
		{
			return new HashSet<String>();
		}
		else
		{
			return spanningFragments;
		}
	}

}
