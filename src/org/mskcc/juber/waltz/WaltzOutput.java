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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * @author Juber Patel
 *
 */
public class WaltzOutput
{
	/**
	 * writers with and without duplicates:
	 * 
	 * 1. pileup
	 * 2. interval/region stats/info
	 * 3. genotypes of interest
	 */

	private String sampleName;
	private BufferedWriter pileupWriter;
	private BufferedWriter pileupWithoutDuplicatesWriter;
	private BufferedWriter intervalsWriter;
	private BufferedWriter intervalsWithoutDuplicatesWriter;
	private BufferedWriter genotypesWriter;
	// private BufferedWriter genotypesWithoutDuplicatesWriter;
	/**
	 * only one writer, (no with or without duplicates)
	 */
	private BufferedWriter signatureIntervalsWriter;

	public WaltzOutput(String sampleName)
	{
		this.sampleName = sampleName;
	}

	public void enableForMetrics() throws IOException
	{
		this.pileupWriter = new BufferedWriter(
				new FileWriter(sampleName + "-pileup.txt"));
		this.pileupWithoutDuplicatesWriter = new BufferedWriter(
				new FileWriter(sampleName + "-pileup-without-duplicates.txt"));
		this.intervalsWriter = new BufferedWriter(
				new FileWriter(sampleName + "-intervals.txt"));
		this.intervalsWithoutDuplicatesWriter = new BufferedWriter(
				new FileWriter(
						sampleName + "-intervals-without-duplicates.txt"));
	}

	public void enableForGenotypes(String mafHeader) throws IOException
	{
		this.genotypesWriter = new BufferedWriter(
				new FileWriter(sampleName + "-genotypes.maf"));
		
		write(genotypesWriter, mafHeader + "\n");
	}

	public void enableForSignatureFinding() throws IOException
	{
		this.signatureIntervalsWriter = new BufferedWriter(
				new FileWriter(sampleName + "-signature-intervals.txt"));
	}

	public void close() throws IOException
	{
		if (pileupWriter != null)
			pileupWriter.close();
		if (pileupWithoutDuplicatesWriter != null)
			pileupWithoutDuplicatesWriter.close();
		if (intervalsWriter != null)
			intervalsWriter.close();
		if (intervalsWithoutDuplicatesWriter != null)
			intervalsWithoutDuplicatesWriter.close();
		if (genotypesWriter != null)
			genotypesWriter.close();
		// if (genotypesWithoutDuplicatesWriter != null)
		// genotypesWithoutDuplicatesWriter.close();
		if (signatureIntervalsWriter != null)
			signatureIntervalsWriter.close();
	}

	public void toPileupWriter(String string) throws IOException
	{
		write(pileupWriter, string);
	}

	public void toPileupWithoutDuplicatesWriter(String string)
			throws IOException
	{
		write(pileupWithoutDuplicatesWriter, string);
	}

	public void toIntervalsWriter(String string) throws IOException
	{
		write(intervalsWriter, string);
	}

	public void toIntervalsWithoutDuplicatesWriter(String string)
			throws IOException
	{
		write(intervalsWithoutDuplicatesWriter, string);
	}

	public void toGenotypesWriter(String string) throws IOException
	{
		write(genotypesWriter, string);
	}

	// public void toGenotypesWithoutDuplicatesWriter(String string)
	// throws IOException
	// {
	// write(genotypesWithoutDuplicatesWriter, string);
	// }

	public void toSignatureIntervalsWriter(String string) throws IOException
	{
		write(signatureIntervalsWriter, string);
	}

	private void write(BufferedWriter writer, String string) throws IOException
	{
		synchronized (writer)
		{
			writer.write(string);
			writer.flush();
		}
	}

	public String getSampleName()
	{
		return sampleName;
	}

}
