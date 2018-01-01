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
package org.mskcc.juber.waltz.pileup.processors;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.mskcc.juber.genotype.Genotype;
import org.mskcc.juber.genotype.GenotypeEventType;
import org.mskcc.juber.genotype.GenotypeID;
import org.mskcc.juber.waltz.WaltzOutput;
import org.mskcc.juber.waltz.pileup.ContineousSpan;
import org.mskcc.juber.waltz.pileup.FragmentSpan;
import org.mskcc.juber.waltz.pileup.RegionPileupView;

import com.google.common.collect.Sets;

public class GenotypingProcessorOld implements PileupProcessor
{
	private RegionPileupView pileup;
	private PileupMetricsProcessor metricsProcessor;
	private List<GenotypeIDWithName> genotypeIDsWithName;

	public GenotypingProcessorOld(File lociFile) throws IOException
	{
		genotypeIDsWithName = new ArrayList<GenotypeIDWithName>();

		BufferedReader reader = new BufferedReader(new FileReader(lociFile));
		String line = null;
		String[] parts = null;

		while ((line = reader.readLine()) != null)
		{
			parts = line.split("\t");
			GenotypeID genotypeID = makeGenotypeID(parts);
			genotypeIDsWithName
					.add(new GenotypeIDWithName(genotypeID, parts[4]));
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
		if (diff == 0 && ref.length == 1)
		{
			// SNV
			return new GenotypeID(GenotypeEventType.SNV, parts[0], position,
					ref, alt);
		}
		else if (diff == 0 && ref.length > 1)
		{
			// MNV
			return new GenotypeID(GenotypeEventType.MNV, parts[0], position,
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
		metricsProcessor = new PileupMetricsProcessor();
		metricsProcessor.setRegionPileupView(pileup);
	}

	@Override
	public void processRegion(WaltzOutput output) throws IOException
	{
		metricsProcessor.processRegion(output);
		processGenotypesOld(output);

		// free the memory once we are done
		this.pileup = null;
	}

	private void processGenotypesOld(WaltzOutput output) throws IOException
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

			// this genotype is not in this region
			if (!pileup.contains(genotypeID))
			{
				continue;
			}

			String outString = null;
			if (genotypeID.type == GenotypeEventType.SNV)
			{
				// process SNV
				outString = processSNV(genotypeID, genotypeName);
			}
			else if (genotypeID.type == GenotypeEventType.MNV)
			{
				// process SNV
				outString = processMNV(genotypeID, genotypeName);
			}
			else if (genotypeID.type == GenotypeEventType.DELETION
					&& genotypeID.ref.length == 2)
			{
				outString = processSinglebaseDeletion(genotypeID, genotypeName);
			}
			else
			{
				// process genotype
				outString = processGenotype(genotypeID, genotypeName);
			}

			output.toGenotypesWriter(outString + "\n");
		}
	}

	private void processGenotypes(WaltzOutput output) throws IOException
	{
		Set<GenotypeIDWithName> regionGenotypes = new HashSet<GenotypeIDWithName>();

		// pick out the genotypes contained in the current region/pileup
		for (int i = 0; i < genotypeIDsWithName.size(); i++)
		{
			GenotypeIDWithName genotypeIDWithName = genotypeIDsWithName.get(i);
			if (genotypeIDWithName == null)
			{
				continue;
			}

			if (pileup.contains(genotypeIDWithName.genotypeID))
			{
				regionGenotypes.add(genotypeIDWithName);
			}
		}

		String outString = null;
		outString = processGenotypeSet(regionGenotypes);
		output.toGenotypesWriter(outString + "\n");
	}

	private String processGenotypeSet(Set<GenotypeIDWithName> regionGenotypes)
	{
		StringBuilder string = null;

		Map<Set<GenotypeIDWithName>, Set<String>> powerSetFragments = getSupportingFragmentsForPowerSet(
				regionGenotypes);

		// process each element of the powerset
		for (Set<GenotypeIDWithName> s : powerSetFragments.keySet())
		{
			Genotype genotype = null;

			// a single, traditional genotype
			if (s.size() == 1)
			{
				Iterator<GenotypeIDWithName> it = s.iterator();
				GenotypeIDWithName genotypeIDWithName = it.next();
				genotype = new Genotype(genotypeIDWithName.name,
						genotypeIDWithName.genotypeID);
			}
			else
			{
				StringBuilder name = null;
				// build name
				for (GenotypeIDWithName genotypeIDWithName : s)
				{
					if (name == null)
					{
						name = new StringBuilder(genotypeIDWithName.name);
					}
					else
					{
						name.append(genotypeIDWithName.name);
					}
				}

				genotype = new Genotype(name.toString());
			}

			ContineousSpan span = computeSpan(s);

			// record spanning fragments
			for (FragmentSpan fragmentSpan : pileup.fragmentSpans.values())
			{
				if (fragmentSpan.spans(span.contig, span.start, span.end))
				{
					genotype.totalCoverage++;

					if (!fragmentSpan.isDuplicate())
					{
						genotype.uniqueCoverage++;
					}
				}
			}

			// ignore genotype sets that have 0 spanning fragments
			if (genotype.totalCoverage == 0)
			{
				continue;
			}

			Set<String> fragments = powerSetFragments.get(s);

			// record supporting fragments
			if (fragments != null)
			{
				genotype.totalSupportingCoverage = fragments.size();
				for (String fragment : fragments)
				{
					if (!pileup.fragmentSpans.get(fragment).isDuplicate())
					{
						genotype.uniqueSupportingCoverage++;
					}
				}
			}

			// record genotype
			if (string == null)
			{
				string = new StringBuilder(genotype.toString());
			}
			else
			{
				string.append(System.lineSeparator());
				string.append(genotype.toString());
			}
		}

		return string.toString();

	}

	/**
	 * compute the total span of the given set of genotype IDs
	 * 
	 * @param s
	 * @return
	 */
	private ContineousSpan computeSpan(Set<GenotypeIDWithName> s)
	{
		String contig = null;
		int start = 0;
		int end = 0;

		for (GenotypeIDWithName genotypeIDWithName : s)
		{
			GenotypeID genotypeID = genotypeIDWithName.genotypeID;

			// first member
			if (contig == null)
			{
				contig = genotypeID.contig;
				start = genotypeID.position;
				end = genotypeID.endPosition;
			}

			// contig check
			if (!contig.equals(genotypeID.contig))
			{
				System.err
						.println("Contigs in a mutation group don't match up!");
				System.err.println(contig + " AND " + genotypeID.toString());
				System.exit(1);
			}

			// see if the span needs to be increased
			if (genotypeID.position < start)
			{
				start = genotypeID.position;
			}

			if (genotypeID.endPosition > end)
			{
				end = genotypeID.endPosition;
			}
		}

		return new ContineousSpan(contig, start, end);

	}

	private String processSNV(GenotypeID genotypeID, String genotypeName)
	{
		int pileupIndex = pileup.getIndex(genotypeID.contig,
				genotypeID.position);
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

		return genotype.toString();
	}

	private String processSinglebaseDeletion(GenotypeID genotypeID,
			String genotypeName)
	{
		int pileupIndex = pileup.getIndex(genotypeID.contig,
				genotypeID.position);
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

		return genotype.toString();
	}

	private String processGenotype(GenotypeID genotypeID, String genotypeName)
	{
		Genotype genotype = new Genotype(genotypeName, genotypeID);

		// find the fragments supporting genotype
		Set<String> fragments = pileup.genotypes.get(genotypeID);

		// record supporting fragments
		if (fragments != null)
		{
			genotype.totalSupportingCoverage = fragments.size();
			for (String fragment : fragments)
			{
				if (!pileup.fragmentSpans.get(fragment).isDuplicate())
				{
					genotype.uniqueSupportingCoverage++;
				}
			}
		}

		// record spanning fragments
		for (FragmentSpan fragmentSpan : pileup.fragmentSpans.values())
		{
			if (fragmentSpan.spans(genotypeID.contig, genotypeID.position,
					genotypeID.endPosition))
			{
				genotype.totalCoverage++;

				if (!fragmentSpan.isDuplicate())
				{
					genotype.uniqueCoverage++;
				}
			}
		}

		/*
		 * // populate total and unique coverage values by taking average
		 * int pileupIndex = pileup.getIndex(genotypeID.contig,
		 * genotypeID.position + 1);
		 * int deletedBases = genotypeID.ref.length;
		 * for (int i = 0; i < deletedBases; i++, pileupIndex++)
		 * {
		 * genotype.totalCoverage += pileup.positions[pileupIndex]
		 * .getCoverage();
		 * genotype.uniqueCoverage +=
		 * pileup.positionsWithoutDuplicates[pileupIndex]
		 * .getCoverage();
		 * }
		 * 
		 * genotype.totalCoverage = genotype.totalCoverage / deletedBases;
		 * genotype.uniqueCoverage = genotype.uniqueCoverage / deletedBases;
		 */

		return genotype.toString();
	}

	private String processMNV(GenotypeID genotypeID, String genotypeName)
	{
		Genotype genotype = new Genotype(genotypeName, genotypeID);

		List<GenotypeID> SNVs = MNVToSNVs(genotypeID);
		Map<Set<GenotypeID>, Set<String>> powerSetFragments = getSupportingFragmentsForPowerSet(
				new HashSet<GenotypeID>(SNVs));
		Set<String> fragments = null;
		// find the full MNV set
		for (Set<GenotypeID> s : powerSetFragments.keySet())
		{
			if (s.size() == genotypeID.ref.length)
			{
				fragments = powerSetFragments.get(s);
			}
		}

		// record supporting fragments
		if (fragments != null)
		{
			genotype.totalSupportingCoverage = fragments.size();
			for (String fragment : fragments)
			{
				if (!pileup.fragmentSpans.get(fragment).isDuplicate())
				{
					genotype.uniqueSupportingCoverage++;
				}
			}
		}

		// record spanning fragments
		for (FragmentSpan fragmentSpan : pileup.fragmentSpans.values())
		{
			if (fragmentSpan.spans(genotypeID.contig, genotypeID.position,
					genotypeID.endPosition))
			{
				genotype.totalCoverage++;

				if (!fragmentSpan.isDuplicate())
				{
					genotype.uniqueCoverage++;
				}
			}
		}

		return genotype.toString();

	}

	/**
	 * Iterate over the powerset of the given set of genotypes and find the
	 * fragments supporting each non-empty element of the powerset
	 * 
	 * @param sNVs
	 * @return
	 */
	private Map<Set<GenotypeIDWithName>, Set<String>> getSupportingFragmentsForPowerSet(
			Set<GenotypeIDWithName> genotypeIDsWithName)
	{
		// find supporting fragments for individual genotypes
		Map<GenotypeIDWithName, Set<String>> supportingFragments = new HashMap<GenotypeIDWithName, Set<String>>();
		for (GenotypeIDWithName genotypeIDWithName : genotypeIDsWithName)
		{
			Set<String> fragments = pileup.genotypes
					.get(genotypeIDWithName.genotypeID);
			supportingFragments.put(genotypeIDWithName, fragments);
		}

		// compute powerset
		Set<Set<GenotypeIDWithName>> powerSet = Sets
				.powerSet(genotypeIDsWithName);

		// populate the returning set
		Map<Set<GenotypeIDWithName>, Set<String>> returningSet = new HashMap<Set<GenotypeIDWithName>, Set<String>>();
		for (Set<GenotypeIDWithName> s : powerSet)
		{
			if (s.isEmpty())
			{
				continue;
			}

			Set<String> fragments = getIntersection(s, supportingFragments);
			returningSet.put(s, fragments);
		}

		return returningSet;
	}

	/**
	 * find the fragments that have all the genotypes in genotypeIDs, using
	 * fragment sets supporting individual genotypeIDs
	 * 
	 * @param genotypeIDs
	 * @param supportingFragments
	 * @return
	 */
	private Set<String> getIntersection(Set<GenotypeIDWithName> genotypeIDs,
			Map<GenotypeIDWithName, Set<String>> supportingFragments)
	{
		Set<String> intersecting = null;
		for (GenotypeIDWithName genotypeIDWithName : genotypeIDs)
		{
			Set<String> supporting = supportingFragments
					.get(genotypeIDWithName);

			if (intersecting == null)
			{
				intersecting = supporting;
			}
			else
			{
				intersecting.retainAll(supporting);
			}
		}

		return intersecting;
	}

	/**
	 * break an MNV into constituent SNVs
	 * 
	 * @param genotypeID
	 * @return
	 */
	private List<GenotypeID> MNVToSNVs(GenotypeID MNVID)
	{
		List<GenotypeID> SNVs = new ArrayList<GenotypeID>();
		for (int i = 0; i < MNVID.ref.length; i++)
		{
			GenotypeID SNV = new GenotypeID(GenotypeEventType.SNV, MNVID.contig,
					MNVID.position + i, new byte[] { MNVID.ref[i] },
					new byte[] { MNVID.alt[i] });
			SNVs.add(SNV);
		}

		return SNVs;
	}

	private String processMultibaseDeletion(GenotypeID genotypeID,
			String genotypeName)
	{
		Genotype genotype = new Genotype(genotypeName, genotypeID);

		// find the fragments supporting genotype
		Set<String> fragments = pileup.genotypes.get(genotypeID);

		// record supporting fragments
		if (fragments != null)
		{
			genotype.totalSupportingCoverage = fragments.size();
			for (String fragment : fragments)
			{
				if (!pileup.fragmentSpans.get(fragment).isDuplicate())
				{
					genotype.uniqueSupportingCoverage++;
				}
			}
		}

		// record spanning fragments
		for (FragmentSpan fragmentSpan : pileup.fragmentSpans.values())
		{
			if (fragmentSpan.spans(genotypeID.contig, genotypeID.position + 1,
					genotypeID.position + genotypeID.alt.length - 1))
			{
				genotype.totalCoverage++;

				if (!fragmentSpan.isDuplicate())
				{
					genotype.uniqueCoverage++;
				}
			}
		}

		/*
		 * // populate total and unique coverage values by taking average
		 * int pileupIndex = pileup.getIndex(genotypeID.contig,
		 * genotypeID.position + 1);
		 * int deletedBases = genotypeID.ref.length;
		 * for (int i = 0; i < deletedBases; i++, pileupIndex++)
		 * {
		 * genotype.totalCoverage += pileup.positions[pileupIndex]
		 * .getCoverage();
		 * genotype.uniqueCoverage +=
		 * pileup.positionsWithoutDuplicates[pileupIndex]
		 * .getCoverage();
		 * }
		 * 
		 * genotype.totalCoverage = genotype.totalCoverage / deletedBases;
		 * genotype.uniqueCoverage = genotype.uniqueCoverage / deletedBases;
		 */

		return genotype.toString();
	}

	private String processInsertion(GenotypeID genotypeID, String genotypeName)
	{
		// find the special genotype
		Genotype genotype = null; // pileup.genotypes.get(genotypeID);
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

		return genotype.toString();
	}

	private class GenotypeIDWithName
	{
		public final GenotypeID genotypeID;
		public final String name;

		public GenotypeIDWithName(GenotypeID genotypeID, String name)
		{
			this.genotypeID = genotypeID;
			this.name = name;
		}
	}

}
