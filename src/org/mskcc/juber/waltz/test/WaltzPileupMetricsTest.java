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
public class WaltzPileupMetricsTest
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

		String module = "PileupMetrics";
		// String module = "Genotyping";
		// String module = "SignatureFinding";

		// String filterType = "BasicFilter";
		String minimumMappingQuality = "1";

		String bamFile = "bamFiles/Pool-cfDNA-30ng-1-5uM-IGO-05500-DL-5_bc219_Pool-05500-DL-Tube3-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";

		// String referenceFasta =
		// "/Users/patelj1/resources/hg19-ucsc/human_hg19.fa";
		String referenceFasta = "/Users/patelj1/resources/hg19-ncbi/Homo_sapiens_assembly19.fasta";

		// TODO re-check if this is the right file!!
		// String bedFile =
		// "/Users/patelj1/workspace/CustomCapture/bedFiles/newBedFiles/BRAF-bad.bed";
		// String intervalsBedFile = "BRAF.bed";
		String intervalsBedFile = "bedFiles/Sarath-10-genes.bed";
		// String intervalsBedFile =
		// "/Users/patelj1/workspace/Marianas/bedFiles/STAG2-CDKN2A-TP53-EWSR1-my-targets.bed";
		// String intervalsBedFile =
		// "/Users/patelj1/workspace/Marianas/bedFiles/ESR1-TP53-PIK3CA-PTEN-ERBB2-AKT1-CDH1-GATA3-my-targets.bed";

		String moduleArgument = null;

		Waltz.main(new String[] { module, minimumMappingQuality, bamFile,
				referenceFasta, intervalsBedFile, moduleArgument });

	}
}