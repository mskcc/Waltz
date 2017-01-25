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
package org.mskcc.juber.waltz;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

/**
 * @author Juber Patel
 * 
 *         static utilities
 *
 */
public class Util
{
	/**
	 * sort properly, on both chr and position
	 * 
	 * @param intervals
	 */
	public static void sort(List<Interval> intervals)
	{
		intervals.sort(new Comparator<Interval>()
		{
			@Override
			public int compare(Interval o1, Interval o2)
			{
				int r1 = o1.getContig().compareTo(o2.getContig());

				if (r1 != 0)
				{
					return r1;
				}
				else
				{
					return o1.getStart() - o2.getStart();
				}
			}
		});
	}

	/**
	 * sort properly, on both chr and position
	 * 
	 * @param intervals
	 */
	public static IntervalList sortAndUniquify(IntervalList intervals)
	{
		List<Interval> list = new ArrayList<Interval>(intervals.size());
		list.addAll(intervals.getIntervals());
		sort(list);
		IntervalList result = new IntervalList(intervals.getHeader());
		result.addall(list);

		// now do unique
		list = IntervalList.getUniqueIntervals(result, false);
		result = new IntervalList(intervals.getHeader());
		result.addall(list);
		return result;
	}
}
