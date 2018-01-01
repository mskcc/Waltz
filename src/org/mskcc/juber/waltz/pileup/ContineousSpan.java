/**
 * 
 */
package org.mskcc.juber.waltz.pileup;

/**
 * @author Juber Patel
 *
 */
public class ContineousSpan
{
	public String contig;
	public int start;
	public int end;

	public ContineousSpan(String contig, int start, int end)
	{
		this.contig = contig;
		this.start = start;
		this.end = end;
	}
}
