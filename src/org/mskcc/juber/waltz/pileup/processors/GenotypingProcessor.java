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
package org.mskcc.juber.waltz.pileup.processors;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.mskcc.juber.genotype.Genotype;
import org.mskcc.juber.genotype.GenotypeEventType;
import org.mskcc.juber.genotype.GenotypeID;
import org.mskcc.juber.waltz.WaltzOutput;
import org.mskcc.juber.waltz.pileup.RegionPileupView;

public class GenotypingProcessor implements PileupProcessor
{
	private RegionPileupView pileup;
	private MetricsProcessor metricsProcessor;
	private List<GenotypeID> genotypeIDs;
	private List<String> genotypeNames;

	public GenotypingProcessor(File lociFile) throws IOException
	{
		genotypeIDs = new ArrayList<GenotypeID>();
		genotypeNames = new ArrayList<String>();

		BufferedReader reader = new BufferedReader(new FileReader(lociFile));
		String line = null;
		String[] parts = null;

		while ((line = reader.readLine()) != null)
		{
			parts = line.split("\t");
			GenotypeID genotypeID = makeGenotypeID(parts);
			genotypeIDs.add(genotypeID);
			genotypeNames.add(parts[4]);
		}

		reader.close();
	}

	/**
	 * make genotype ids from a parsed line from an input mutations files
	 * 
	 * @param parts
	 * @return
	 */
	private GenotypeID makeGenotypeID(String[] parts)
	{
		int position = Integer.parseInt(parts[1]);
		byte[] ref = new byte[parts[2].length()];
		for (int i = 0; i < ref.length; i++)
		{
			ref[i] = (byte) parts[2].charAt(i);
		}

		byte[] alt = new byte[parts[3].length()];
		for (int i = 0; i < alt.length; i++)
		{
			alt[i] = (byte) parts[3].charAt(i);
		}

		int diff = ref.length - alt.length;
		if (diff == 0 && parts[2].length() == 1)
		{
			// SNV
			return new GenotypeID(GenotypeEventType.SNV, parts[0], position,
					ref, alt);
		}
		else if (diff > 0)
		{
			// deletion
			return new GenotypeID(GenotypeEventType.DELETION, parts[0],
					position, ref, alt);
		}
		else if (diff < 0)
		{
			// insertion
			return new GenotypeID(GenotypeEventType.INSERTION, parts[0],
					position, ref, alt);
		}
		else
		{
			System.err.println("Not Supported: " + parts[0] + "\t" + parts[1]
					+ "\t" + parts[2] + "\t" + parts[3] + "\t" + parts[4]
					+ "\n");
			return null;
		}
	}

	@Override
	public void setRegionPileupView(RegionPileupView view)
	{
		this.pileup = view;
		metricsProcessor = new MetricsProcessor();
		metricsProcessor.setRegionPileupView(pileup);
	}

	@Override
	public void processRegion(WaltzOutput output) throws IOException
	{
		metricsProcessor.processRegion(output);
		processGenotypes(output);

		// free the memory once we are done
		this.pileup = null;
	}

	private void processGenotypes(WaltzOutput output) throws IOException
	{
		// process the genotypes list one by one
		for (int i = 0; i < genotypeIDs.size(); i++)
		{
			GenotypeID genotypeID = genotypeIDs.get(i);
			if (genotypeID == null)
			{
				continue;
			}

			String genotypeName = genotypeNames.get(i);
			int pileupIndex = pileup.getIndex(genotypeID.contig,
					genotypeID.position);

			// this genotype is not in this region
			if (pileupIndex == -1)
			{
				continue;
			}

			String outString = null;
			if (genotypeID.type == GenotypeEventType.SNV)
			{
				// process SNV
				outString = processSNV(genotypeID, genotypeName, pileupIndex);
			}
			else if (genotypeID.type == GenotypeEventType.DELETION)
			{
				// process deletion
				if (genotypeID.ref.length == 2)
				{
					outString = processSinglebaseDeletion(genotypeID,
							genotypeName, pileupIndex);
				}
				else
				{
					outString = processMultibaseDeletion(genotypeID,
							genotypeName);
				}
			}
			else if (genotypeID.type == GenotypeEventType.INSERTION)
			{
				// process insertion
				outString = processInsertion(genotypeID, genotypeName);
			}

			output.toGenotypesWriter(outString + "\n");
		}
	}

	private String processSNV(GenotypeID genotypeID, String genotypeName,
			int pileupIndex)
	{
		byte refFromPileup = pileup.positions[pileupIndex].getRefBase();

		if (genotypeID.ref[0] != refFromPileup)
		{
			System.err.println(
					"Ref base for the locus does not match the ref base in the pileup!");
			System.err.println(genotypeID.toString());
			System.err.println("Ref in pileup: " + (char) refFromPileup);
			System.exit(1);
		}

		// create and populate the genotype with correct values
		Genotype genotype = new Genotype(genotypeName, genotypeID);
		genotype.totalCoverage = pileup.positions[pileupIndex].getCoverage();
		genotype.totalSupportingCoverage = pileup.positions[pileupIndex]
				.getCount(genotypeID.alt[0]);
		genotype.uniqueCoverage = pileup.positionsWithoutDuplicates[pileupIndex]
				.getCoverage();
		genotype.uniqueSupportingCoverage = pileup.positionsWithoutDuplicates[pileupIndex]
				.getCount(genotypeID.alt[0]);

		return genotype.toString(genotypeID, genotypeName);
	}

	private String processSinglebaseDeletion(GenotypeID genotypeID,
			String genotypeName, int pileupIndex)
	{
		byte refFromPileup = pileup.positions[pileupIndex].getRefBase();
		byte deletionBase = (byte) 'D';

		if (genotypeID.ref[0] != refFromPileup)
		{
			System.err.println(
					"Ref base for the locus does not match the ref base in the pileup!");
			System.err.println(genotypeID.toString());
			System.err.println("Ref in pileup: " + (char) refFromPileup);
			System.exit(1);
		}

		// move to the actual position of deletion
		pileupIndex++;

		// create and populate the genotype with correct values
		Genotype genotype = new Genotype(genotypeName, genotypeID);
		genotype.totalCoverage = pileup.positions[pileupIndex].getCoverage();
		genotype.totalSupportingCoverage = pileup.positions[pileupIndex]
				.getCount(deletionBase);
		genotype.uniqueCoverage = pileup.positionsWithoutDuplicates[pileupIndex]
				.getCoverage();
		genotype.uniqueSupportingCoverage = pileup.positionsWithoutDuplicates[pileupIndex]
				.getCount(deletionBase);

		return genotype.toString(genotypeID, genotypeName);
	}

	private String processMultibaseDeletion(GenotypeID genotypeID,
			String genotypeName)
	{
		// find the special genotype
		Genotype genotype = pileup.specialGenotypes.get(genotypeID);
		if (genotype == null)
		{
			// this genotype was not recorded, create empty genotype
			genotype = new Genotype(genotypeName, genotypeID);
		}

		// populate total and unique coverage values by taking average
		int pileupIndex = pileup.getIndex(genotypeID.contig,
				genotypeID.position + 1);
		int deletedBases = genotypeID.ref.length;
		for (int i = 0; i < deletedBases; i++, pileupIndex++)
		{
			genotype.totalCoverage += pileup.positions[pileupIndex]
					.getCoverage();
			genotype.uniqueCoverage += pileup.positionsWithoutDuplicates[pileupIndex]
					.getCoverage();
		}

		genotype.totalCoverage = genotype.totalCoverage / deletedBases;
		genotype.uniqueCoverage = genotype.uniqueCoverage / deletedBases;

		return genotype.toString(genotypeID, genotypeName);
	}

	private String processInsertion(GenotypeID genotypeID, String genotypeName)
	{
		// find the special genotype
		Genotype genotype = pileup.specialGenotypes.get(genotypeID);
		if (genotype == null)
		{
			// this genotype was not recorded, create empty genotype
			genotype = new Genotype(genotypeName, genotypeID);
		}

		// populate total and unique coverage values
		// use the coverage at preceding genomic position as proxy
		int pileupIndex = pileup.getIndex(genotypeID.contig,
				genotypeID.position);
		genotype.totalCoverage = pileup.positions[pileupIndex].getCoverage();
		genotype.uniqueCoverage = pileup.positionsWithoutDuplicates[pileupIndex]
				.getCoverage();

		return genotype.toString(genotypeID, genotypeName);
	}

}
