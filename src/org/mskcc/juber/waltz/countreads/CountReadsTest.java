/**
 * 
 */
package org.mskcc.juber.waltz.countreads;

import java.io.IOException;

import org.mskcc.juber.util.CustomCaptureException;

/**
 * @author Juber Patel
 *
 */
public class CountReadsTest
{

	/**
	 * @param args
	 * @throws IOException
	 * @throws CustomCaptureException
	 */
	public static void main(String[] args)
			throws IOException, CustomCaptureException
	{
		String geneList = "/Users/patelj1/resources/gene-list/juber-hg19-gene-list.bed";
		String coverageThreshold = "10";
		
		
		//String bam = "/Users/patelj1/workspace/PUMA/5500-AR/FinalBams/DS-puma-0027-PL-C3-IGO-05500-AR-33_bc67_5500-AR_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		//String bedFile = "bedFiles/ERBB2.bed";
		
		//String bam = "bamFiles/MCC_P-0014336-T01_IGO_05500_DG_11_S79_L004.bam";
		// String bam = "bamFiles/Pool-cfDNA-30ng-1-5uM-IGO-05500-DL-5_bc219_Pool-05500-DL-Tube3-1_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		//String bam = "bamFiles/Pan-Cancer-F1_S14_L001.bam";
		//String bam = "t.bam";
		//String bam = "bamFiles/DS-puma-0006-PL-C17-IGO-05500-BO-7_bc107_5500-BO_L000_mrg_cl_aln_srt_MD_IR_FX_BR.bam";
		String bam = "bamFiles/"
				+ "C-WUC3WV-N002-d_cl_aln_srt_MD_IR_FX_BR.bam";

		// String bedFile = "bedFiles/impact410-mcpyv-ebv-hpv.bed";
		// String bedFile = "bedFiles/Sarath-10-genes.bed";
		String bedFile = "bedFiles/pan-cancer-panel.bed";
		//String bedFile = "bedFiles/ERBB2-PUMA.bed";

		CountReads.main(
				new String[] { bam, coverageThreshold, geneList, bedFile });
	}
}
