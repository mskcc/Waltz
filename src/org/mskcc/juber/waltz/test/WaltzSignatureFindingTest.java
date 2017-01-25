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
package org.mskcc.juber.waltz.test;

import java.io.IOException;
import java.util.concurrent.ExecutionException;

import org.mskcc.juber.waltz.Waltz;

/**
 * @author Juber Patel
 *
 */
public class WaltzSignatureFindingTest
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
		// String module = "Genotyping";
		String module = "SignatureFinding";

		// String filterType = "BasicFilter";
		String filterType = "QualityFilter";

		// String bamFile =
		// "/Users/patelj1/workspace/Shukla/FinalBams/ES-CTDNA-15-01-IGO-05500-AQ-4_bc42_5500-AQ_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		// String bamFile =
		// "/Users/patelj1/workspace/Shukla/FinalBams/ES-CTDNA-09-01-IGO-05500-AQ-17_bc79_5500-AQ_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		String bamFile = "/Users/patelj1/workspace/Shukla/run-5500-AQ/FinalBams/ES-CTDNA-03-01-IGO-05500-AQ-9_bc49_5500-AQ_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		// String bamFile =
		// "/Users/patelj1/workspace/Moynahan/FinalBams/1196-2-IGO-05500-AL-21_bc37_5500-AL_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";

		// String referenceFasta =
		// "/Users/patelj1/resources/hg19-ucsc/human_hg19.fa";
		String referenceFasta = "/Users/patelj1/resources/hg19-ncbi/Homo_sapiens_assembly19.fasta";

		// TODO re-check if this is the right file!!
		// String bedFile =
		// "/Users/patelj1/workspace/CustomCapture/bedFiles/newBedFiles/BRAF-bad.bed";
		// String intervalsBedFile = "BRAF.bed";
		String intervalsBedFile = "/Users/patelj1/workspace/Marianas/bedFiles/EWSR1.bed";
		// String intervalsBedFile =
		// "/Users/patelj1/workspace/Marianas/bedFiles/STAG2-CDKN2A-TP53-EWSR1.bed";
		// String intervalsBedFile =
		// "/Users/patelj1/workspace/Marianas/bedFiles/ESR1-TP53-PIK3CA-PTEN-ERBB2-AKT1-CDH1-GATA3.bed";
		// String lociBedFile = "BRAFLoci.bed";

		String signaturesToLookFor = "TranslocationBreakpoint";

		Waltz.main(new String[] { module, filterType, bamFile,
				referenceFasta, intervalsBedFile, signaturesToLookFor });
	}
}
