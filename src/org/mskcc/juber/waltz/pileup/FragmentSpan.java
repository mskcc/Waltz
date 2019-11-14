/**
 * 
 */
package org.mskcc.juber.waltz.pileup;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMRecord;

/**
 * @author Juber Patel
 * 
 *         Span of a fragment, including read1 and read2 in case of paired end
 *         sequencing. It could be more than 2 records. But records are assumed
 *         to be sorted by alignment start, as in a position sorted bam file.
 *         Each record will be added independently since mutations from CIGAR
 *         strings cannot span records.
 *
 */
public class FragmentSpan
{
	private List<ContineousSpan> recordSpans;
	private boolean duplicate;

	public FragmentSpan(boolean duplicate)
	{
		this.duplicate = duplicate;
		recordSpans = new ArrayList<ContineousSpan>();
	}

	public boolean isDuplicate()
	{
		return duplicate;
	}

	/**
	 * add the span and return adjusted, non-overlapping span
	 * 
	 * @param contig
	 * @param start
	 * @param end
	 * @return
	 */
	public ContineousSpan addAndAdjust(SAMRecord record)
	{
		String contig = record.getContig();
		int start = record.getAlignmentStart();
		int end = record.getAlignmentEnd();
		ContineousSpan newRecordSpan = new ContineousSpan(contig, start, end);

		// Calculate the adjusted, non-overlapping part of the given record
		// Since intervals will be fed in position sorted manner, there are only
		// 4 possibilities vis a vis a recordSpan:
		// 1. given interval is not on the same chromosome as recordSpan
		// 2. given interval is completely contained in the recordSpan
		// 3. given interval is to the right of and completely disjoint from the
		// recordSpan.
		// 4. only a left part of the given interval overlaps with the
		// recordSpan.

		for (ContineousSpan recordSpan : recordSpans)
		{
			// not on same contig, don't bother
			if (!contig.equals(recordSpan.contig))
			{
				continue;
			}

			// fully contained in current record span
			else if (start >= recordSpan.start && end <= recordSpan.end)
			{
				// we are done
				start = end + 1;
				break;
			}

			// to the right and disjoint from the current record span
			else if (start > recordSpan.end)
			{
				continue;
			}

			// left side of the interval is partially contained in the current
			// record span
			// clip the overlapping part
			else if (recordSpan.end >= start && recordSpan.end <= end)
			{
				start = recordSpan.end + 1;

			}
		}

		// add the given record's span
		recordSpans.add(newRecordSpan);

		// return adjusted span
		if (start <= end)
		{
			return new ContineousSpan(contig, start, end);
		}

		return null;
	}

	/**
	 * Does this fragment completely span the given interval? Must be fully
	 * spanned by at least one record.
	 * 
	 * @param contig
	 * @param start
	 * @param end
	 * @return
	 */
	public int spanningReads(String contig, int start, int end)
	{
		int spanningReads = 0;

		for (ContineousSpan recordSpan : recordSpans)
		{
			// not on same contig, don't bother
			if (!contig.equals(recordSpan.contig))
			{
				continue;
			}

			// fully contained in current record span
			if (start >= recordSpan.start && end <= recordSpan.end)
			{
				spanningReads++;
			}
		}

		return spanningReads;
	}
}
