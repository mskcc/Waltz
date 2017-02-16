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
package org.mskcc.juber.waltz.pileup.processors.signatures;

import java.util.ArrayList;
import java.util.List;

import org.mskcc.juber.waltz.pileup.PositionPileup;
import org.mskcc.juber.waltz.pileup.RegionPileupView;

/**
 * @author Juber Patel
 * 
 *         Genomic signature of a translocation breakpoint
 *
 */
public class TranslocationBreakpointSignature implements PileupSignature
{

	private static final double threshold = 0.01;
	private static final int minTotalReads = 100;

	@Override
	public List<SignatureLocus> findIn(RegionPileupView pileup)
	{
		List<SignatureLocus> loci = new ArrayList<SignatureLocus>();

		PositionPileup[] positions = pileup.positionsWithoutDuplicates;
		int lastValidPositionIndex = pileup.interval.getEnd()
				- pileup.interval.getStart();
		for (int i = 0; i <= lastValidPositionIndex; i++)
		{
			int total = positions[i].getCoverage();
			if (total < minTotalReads)
			{
				continue;
			}

			String description = null;
			int clips = 0;
			if (positions[i].getHardClipStarts() > positions[i]
					.getHardClipEnds())
			{
				clips = positions[i].getHardClipStarts();
				description = "3' Translocation Breakpoint";
			}
			else
			{
				clips = positions[i].getHardClipEnds();
				description = "5' Translocation Breakpoint";
			}

			// for now, just the ratio of hard-clip-start supporting reads to
			// total reads at the position right before the hard clip (which is
			// called the hard clip start position in the pileup)
			double evidence = (clips * 1.0) / total;

			if (evidence < threshold)
			{
				continue;
			}

			int start = pileup.interval.getStart() + i;
			loci.add(new SignatureLocus(start, start, description,
					(evidence + "\t" + (clips + "/" + total))));
		}

		return loci;
	}

	@Override
	public String getName()
	{
		return "TranslocationBreakpoint";
	}

}
