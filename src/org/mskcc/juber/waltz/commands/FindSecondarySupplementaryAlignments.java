/**
 * 
 */
package org.mskcc.juber.waltz.commands;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * @author Juber Patel
 * 
 *         change all the base qualities in all the alignments in a bam to a
 *         given constant
 *
 */
public class FindSecondarySupplementaryAlignments
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		File inBam = new File(
				"bamFiles/marianas-collapsed/MSK-L-017-cf-IGO-05500-DY-18_bc209_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam");

		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(inBam);
		SAMRecordIterator iterator = reader.iterator();

		long totalMapped = 0;
		long notPrimary = 0;
		long supplementary = 0;

		while (iterator.hasNext())
		{
			SAMRecord record = iterator.next();

			if (record.getReadUnmappedFlag())
			{
				continue;
			}

			totalMapped++;

			if (record.getNotPrimaryAlignmentFlag())
			{
				notPrimary++;
			}

			if (record.getSupplementaryAlignmentFlag())
			{
				supplementary++;
			}
		}

		iterator.close();
		reader.close();

		System.out.println("Mapped: " + totalMapped);
		System.out.println("Not Primary: " + notPrimary + " ("
				+ (notPrimary * 1.0) / totalMapped + ")");
		System.out.println("Supplementary: " + supplementary + " ("
				+ (supplementary * 1.0) / totalMapped + ")");

	}

}
