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
package org.mskcc.juber.waltz.pileup.processors;

import java.io.IOException;

import org.mskcc.juber.waltz.WaltzOutput;
import org.mskcc.juber.waltz.pileup.Fragment;
import org.mskcc.juber.waltz.pileup.PositionPileup;
import org.mskcc.juber.waltz.pileup.RegionPileupView;

/**
 * @author Juber Patel
 *
 */
public class PileupMetricsProcessor implements PileupProcessor
{
	private RegionPileupView pileup;

	@Override
	public void setRegionPileupView(RegionPileupView view)
	{
		this.pileup = view;
	}

	@Override
	public void processRegion(WaltzOutput output) throws IOException
	{
		processPileup(output);
		processInterval(output);

		// free the memory once we are done
		this.pileup = null;
	}

	/**
	 * process the pileup-level info
	 * 
	 * @param outList
	 * @param pileupWriter
	 * @throws IOException
	 */
	private void processPileup(WaltzOutput output) throws IOException
	{
		// pileup that has never been prepped and used
		if (pileup.interval == null)
		{
			return;
		}

		StringBuilder pileupOut = new StringBuilder();
		StringBuilder pileupWithoutDuplicatesOut = new StringBuilder();

		// iterate over positions and collect numbers
		for (int i = 0; i <= pileup.lastValidPositionIndex; i++)
		{
			addToPileup(pileup.positions, i, pileupOut);
			addToPileup(pileup.positionsWithoutDuplicates, i,
					pileupWithoutDuplicatesOut);
		}

		// write pileups to the pileup files
		output.toPileupWriter(pileupOut.toString());
		output.toPileupWithoutDuplicatesWriter(
				pileupWithoutDuplicatesOut.toString());
	}

	private void addToPileup(PositionPileup[] positions, int index,
			StringBuilder pileupOut)
	{
		int pos = pileup.interval.getStart() + index;
		pileupOut.append(pileup.interval.getContig()).append("\t").append(pos)
				.append("\t").append(positions[index].toString()).append("\n");
	}

	/**
	 * add interval level info for the interval of the set pileup
	 * 
	 * @param outList
	 * @param pileupWriter
	 * @throws IOException
	 */
	private void processInterval(WaltzOutput output) throws IOException
	{
		// pileup that has never been prepped and used
		if (pileup.interval == null)
		{
			return;
		}

		// collect and write interval-level metrics
		IntervalMetrics intervalMetrics = new IntervalMetrics(pileup.interval,
				pileup.positions, pileup.fragmentSpans.size());
		output.toIntervalsWriter(intervalMetrics.toString());

		// count unique fragments
		int uniqueFragments = 0;
		for (Fragment fragment : pileup.fragmentSpans.values())
		{
			if (!fragment.isDuplicate())
			{
				uniqueFragments++;
			}
		}

		IntervalMetrics intervalMetricsWithoutDuplicates = new IntervalMetrics(
				pileup.interval, pileup.positionsWithoutDuplicates,
				uniqueFragments);
		output.toIntervalsWithoutDuplicatesWriter(
				intervalMetricsWithoutDuplicates.toString());
	}

}
