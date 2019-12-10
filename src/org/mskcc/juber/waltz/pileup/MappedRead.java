/**
 * 
 */
package org.mskcc.juber.waltz.pileup;

import java.util.HashSet;
import java.util.Set;

import org.mskcc.juber.genotype.GenotypeID;

/**
 * @author Juber Patel
 *
 */
public class MappedRead
{
	public final String contig;
	public final int start;
	public final int end;
	private Set<GenotypeID> genotypes;

	public MappedRead(String contig, int start, int end)
	{
		this.contig = contig;
		this.start = start;
		this.end = end;
		this.genotypes = new HashSet<GenotypeID>();

	}

	public void addGenotype(GenotypeID genotypeID)
	{
		genotypes.add(genotypeID);
	}

	public Set<GenotypeID> getGenotypes()
	{
		return genotypes;
	}
}
