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
package org.mskcc.juber.waltz.pileup.processors;

import org.mskcc.juber.waltz.pileup.PositionPileup;

import htsjdk.samtools.util.Interval;

/**
 * @author Juber Patel
 *
 */
public class IntervalMetrics
{
	private Interval interval;
	private int length;
	private double gcContent;
	private double averageCoverage;
	private int peakCoverage;

	public IntervalMetrics(Interval interval, PositionPileup[] positions)
	{
		this.interval = interval;
		this.length = interval.getEnd() - interval.getStart() + 1;
		int lastValidPositionIndex = interval.getEnd() - interval.getStart();

		// iterate over position pileups
		for (int i = 0; i <= lastValidPositionIndex; i++)
		{
			PositionPileup positionPileup = positions[i];

			byte refBase = positionPileup.getRefBase();
			if (refBase == 'C' || refBase == 'c' || refBase == 'G'
					|| refBase == 'g')
			{
				gcContent = gcContent + 1;
			}

			int coverage = positionPileup.getCoverage();
			averageCoverage = averageCoverage + coverage;

			if (peakCoverage < coverage)
			{
				peakCoverage = coverage;
			}
		}

		// divide
		gcContent = gcContent / length;
		averageCoverage = averageCoverage / length;
	}

	public String toString()
	{
		StringBuilder intervalOutput = new StringBuilder();
		intervalOutput.append(interval.getContig()).append("\t")
				.append(interval.getStart()).append("\t")
				.append(interval.getEnd()).append("\t")
				.append(interval.getName()).append("\t").append(length)
				.append("\t").append(peakCoverage).append("\t")
				.append(averageCoverage).append("\t").append(gcContent)
				.append("\n");
		return intervalOutput.toString();
	}
}
