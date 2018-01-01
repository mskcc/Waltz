/**
 * 
 */
package org.mskcc.juber.waltz.commands;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
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
public class ChangeBamBaseQualities
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		byte baseQuality = 12;
		char BIQuality = '-';

		File inBam = new File(
				"bamFiles/MSK-L-017-cf-IGO-05500-DY-18_bc209_5500-DY-4_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam");
		File outBam = new File(
				"bamFiles/MSK-L-017-cf-bq" + baseQuality + ".bam");

		SamReaderFactory factory = SamReaderFactory.makeDefault();
		SamReader reader = factory.open(inBam);
		SAMFileHeader header = reader.getFileHeader();
		SAMRecordIterator iterator = reader.iterator();

		SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header,
				true, outBam);

		char[] chars = new char[500];
		for (int i = 0; i < chars.length; i++)
		{
			chars[i] = BIQuality;
		}

		byte[] quals = new byte[500];
		for (int i = 0; i < quals.length; i++)
		{
			quals[i] = baseQuality;
		}

		while (iterator.hasNext())
		{
			SAMRecord record = iterator.next();

			// change base qualities, BD, BI, BQ, OQ

			record.setBaseQualities(
					Arrays.copyOf(quals, record.getBaseQualities().length));

			int tagValueLength = record.getStringAttribute("BD").length();
			record.setAttribute("BD", new String(chars, 0, tagValueLength));

			String attribute = record.getStringAttribute("BI");
			if (attribute != null)
			{
				tagValueLength = attribute.length();
				record.setAttribute("BI", new String(chars, 0, tagValueLength));
			}

			attribute = record.getStringAttribute("BQ");
			if (attribute != null)
			{
				tagValueLength = attribute.length();
				record.setAttribute("BQ", new String(chars, 0, tagValueLength));
			}

			attribute = record.getStringAttribute("OQ");
			if (attribute != null)
			{
				tagValueLength = attribute.length();
				record.setAttribute("OQ", new String(chars, 0, tagValueLength));
			}

			writer.addAlignment(record);
		}

		writer.close();
		iterator.close();
		reader.close();

	}

}
