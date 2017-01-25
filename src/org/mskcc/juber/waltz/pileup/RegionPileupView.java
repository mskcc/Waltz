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

import java.util.Map;

import org.mskcc.juber.genotype.Genotype;
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
	public final Map<GenotypeID, Genotype> specialGenotypes;
	public final int insertMin;
	public final int insertMax;
	
	public RegionPileupView(byte[] referenceBases, Interval interval,
			int lastValidPositionIndex, PositionPileup[] positions,
			PositionPileup[] positionsWithoutDuplicates, Map<GenotypeID, Genotype> specialGenotypes, int insertMin,
			int insertMax)
	{
		this.referenceBases = referenceBases;
		this.interval = interval;
		this.lastValidPositionIndex = lastValidPositionIndex;
		this.positions = positions;
		this.positionsWithoutDuplicates = positionsWithoutDuplicates;
		this.specialGenotypes = specialGenotypes;
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

}
