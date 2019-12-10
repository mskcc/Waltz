/**
 * 
 */
package org.mskcc.juber.waltz.pileup;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.mskcc.juber.genotype.GenotypeID;

import htsjdk.samtools.SAMRecord;

/**
 * @author Juber Patel
 * 
 *         A mapped fragment, including read1 and read2 in case of paired end
 *         sequencing. It could be more than 2 records. But records are assumed
 *         to be sorted by alignment start, as in a position sorted bam file.
 *         Each record will be added independently since mutations from CIGAR
 *         strings cannot span records.
 *
 */
public class Fragment
{
	private List<MappedRead> reads;
	private boolean duplicate;

	public Fragment(boolean duplicate)
	{
		this.duplicate = duplicate;
		reads = new ArrayList<MappedRead>();
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
	public MappedRead add(SAMRecord record)
	{
		MappedRead read = new MappedRead(record.getContig(),
				record.getAlignmentStart(), record.getAlignmentEnd());

		// add the new read
		reads.add(read);

		return read;
	}

	/**
	 * return the set of genotypes found in this fragment after applying
	 * readPairMismatchPolicy
	 * 
	 * @param readPairMismtachPolicy
	 * @return
	 */
	public Set<GenotypeID> getGenotypes(int readPairMismtachPolicy)
	{
		Set<GenotypeID> genotypes = new HashSet<GenotypeID>();

		// add all genotypes from all reads
		for (MappedRead read : reads)
		{
			for (GenotypeID genotypeID : read.getGenotypes())
			{
				genotypes.add(genotypeID);
			}
		}

		// remove the ones not allowed by readPairMismatchPolicy
		Iterator<GenotypeID> it = genotypes.iterator();
		while (it.hasNext())
		{
			GenotypeID genotypeID = it.next();

			// remove any genotype with mismatch between reads
			// but keep the genotypes with N
			if (readPairMismtachPolicy == 0 && !altHasN(genotypeID)
					&& rejectedByARead(genotypeID))
			{
				it.remove();
			}
			else if (readPairMismtachPolicy == 1 && altHasN(genotypeID))
			{
				// remove any genotype with N in it
				it.remove();
			}
		}

		return genotypes;
	}

	/**
	 * figure out if there is a read that covers the genotype even partially but
	 * rejects it
	 * 
	 * @param genotypeID
	 * @return
	 */
	private boolean rejectedByARead(GenotypeID genotypeID)
	{
		for (MappedRead read : reads)
		{
			Set<GenotypeID> readGenotypes = read.getGenotypes();

			// genotype fully contained in read.
			if (read.contig.equals(genotypeID.contig)
					&& read.start <= genotypeID.position
					&& read.end >= genotypeID.endPosition)
			{
				// should be there but is not there
				if (!readGenotypes.contains(genotypeID))
				{
					return true;
				}
			}
			else
			{
				// get partial genotype
				GenotypeID partial = genotypeID.partialGenotype(read.contig,
						read.start, read.end);

				// no overlap
				if (partial == null)
				{
					continue;
				}
				else if (!readGenotypes.contains(partial))
				{
					return true;
				}
			}
		}

		return false;
	}

	private boolean altHasN(GenotypeID genotypeID)
	{
		for (int i = 0; i < genotypeID.alt.length; i++)
		{
			if ((char) genotypeID.alt[i] == 'N')
			{
				return true;
			}
		}

		return false;
	}

	/**
	 * 
	 * Calculate the adjusted, non-overlapping part of the given record
	 * Since intervals will be fed in position sorted manner, there are only
	 * 4 possibilities vis a vis a mappedRead:
	 * 1. given interval is not on the same chromosome as recordSpan
	 * 2. given interval is completely contained in the recordSpan
	 * 3. given interval is to the right of and completely disjoint from the
	 * recordSpan.
	 * 4. only a left part of the given interval overlaps with the
	 * recordSpan.
	 * 
	 * Not used anymore!
	 * 
	 * 
	 * @param newRecordSpan
	 * @return
	 */
	private MappedRead getNonOverlappingPartOld(String contig, int start,
			int end)
	{
		for (MappedRead read : reads)
		{
			// not on same contig, don't bother
			if (!contig.equals(read.contig))
			{
				continue;
			}

			// fully contained in current record span
			else if (start >= read.start && end <= read.end)
			{
				// we are done
				start = end + 1;
				break;
			}

			// to the right and disjoint from the current record span
			else if (start > read.end)
			{
				continue;
			}

			// left side of the interval is partially contained in the current
			// record span
			// clip the overlapping part
			else if (read.end >= start && read.end <= end)
			{
				start = read.end + 1;

			}
		}

		return new MappedRead(contig, start, end);

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
	private int spanningReads(String contig, int start, int end)
	{
		int spanningReads = 0;

		for (MappedRead recordSpan : reads)
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
