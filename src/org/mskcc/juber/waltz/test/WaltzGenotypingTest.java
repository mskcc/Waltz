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
package org.mskcc.juber.waltz.test;

import java.io.IOException;
import java.util.concurrent.ExecutionException;

import org.mskcc.juber.waltz.Waltz;

/**
 * @author Juber Patel
 *
 */
public class WaltzGenotypingTest
{

	/**
	 * @param args
	 * @throws ExecutionException
	 * @throws InterruptedException
	 * @throws IOException
	 */
	public static void main(String[] args)
			throws IOException, InterruptedException, ExecutionException
	{
		// TODO a bit more accuracy needed. We are probably missing some reads.

		// String module = "Metrics";
		String module = "Genotyping";

		// String module = "SignatureFinding";

		// String filterType = "BasicFilter";
		String minimumMappingQuality = "20";

		String bamFile = "/Users/patelj1/workspace/Waltz/bamFiles/Hapmap-ctrl-1_S2_001_cl_aln_srt_MD_IR_FX_BR.bam";
		// String bamFile =
		// "/Volumes/innovation/Innovation/projects/Juber/HiSeq/5500-DY/standard-bam-links/"
		// +
		// "PR-pt09-CF01-IGO-05500-DY-27_bc224_5500-DY-5_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";

		String referenceFasta = "/Users/patelj1/resources/hg19-ncbi/Homo_sapiens_assembly19.fasta";

		String intervalsBedFile = "/Users/patelj1/workspace/Waltz/bedFiles/MSK-ACCESS-V1_0-A.bed";
		// String intervalsBedFile = "bedFiles/multibase-substitution-test.bed";

		String moduleArgument = "/Users/patelj1/workspace/Waltz/Genotyping/mutationsFiles/test.maf";
		// String moduleArgument =
		// "/Users/patelj1/workspace/Waltz/Genotyping/mutationsFiles/breast-carcinoma-8-gene-hotspots.txt";

		Waltz.main(new String[] { module, minimumMappingQuality, bamFile,
				referenceFasta, intervalsBedFile, moduleArgument });
	}
}
