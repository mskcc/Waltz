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
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.mskcc.juber.genotype.Genotype;
import org.mskcc.juber.genotype.GenotypeEventType;
import org.mskcc.juber.genotype.GenotypeID;
import org.mskcc.juber.waltz.WaltzOutput;
import org.mskcc.juber.waltz.pileup.ContineousSpan;
import org.mskcc.juber.waltz.pileup.FragmentSpan;
import org.mskcc.juber.waltz.pileup.RegionPileupView;

import com.google.common.collect.Sets;

public class GenotypingProcessor implements PileupProcessor
{
	private RegionPileupView pileup;
	private PileupMetricsProcessor metricsProcessor;
	private List<GenotypeIDWithName> genotypeIDsWithName;
	private Comparator<Entry<GenotypeIDWithName, Set<String>>> fragmentCountComparator;

	public GenotypingProcessor(File lociFile) throws IOException
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

		// create fragment count comparator
		fragmentCountComparator = new Comparator<Entry<GenotypeIDWithName, Set<String>>>()
		{

			@Override
			public int compare(Entry<GenotypeIDWithName, Set<String>> arg0,
					Entry<GenotypeIDWithName, Set<String>> arg1)
			{
				return arg1.getValue().size() - arg0.getValue().size();
			}
		};
	}

	public List<GenotypeID> getGenotypeIDs()
	{
		List<GenotypeID> list = new ArrayList<GenotypeID>();
		for (GenotypeIDWithName g : genotypeIDsWithName)
		{
			list.add(g.genotypeID);
		}

		return list;

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
		// metricsProcessor = new PileupMetricsProcessor();
		// metricsProcessor.setRegionPileupView(pileup);
	}

	@Override
	public void processRegion(WaltzOutput output) throws IOException
	{
		// metricsProcessor.processRegion(output);
		processGenotypes(output);

		// free the memory once we are done
		this.pileup = null;
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

		// no genotypes to be profiled in the current region
		if (regionGenotypes.isEmpty())
		{
			return;
		}

		String outString = null;
		outString = processGenotypeSet(regionGenotypes);
		output.toGenotypesWriter(outString + "\n");
	}

	/**
	 * profile the given genotype set that has all the genotypes for the current
	 * region. Since the powerset of this set can explode quickly, for the time
	 * being we will include in the output:
	 * 
	 * 1. all individual genotypes
	 * 2. all combinations (power set) of genotypes with non-zero supporting
	 * fragments
	 * 
	 * If non-zero genotypes are > 5, then top 5 by number of supporting
	 * fragments are chosen to take the power set
	 * 
	 * 
	 * @param regionGenotypes
	 * @return
	 */
	private String processGenotypeSet(Set<GenotypeIDWithName> regionGenotypes)
	{
		StringBuilder result = null;

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
						name.append("-");
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
			if (result == null)
			{
				result = new StringBuilder(genotype.toString());
			}
			else
			{
				result.append(System.lineSeparator());
				result.append(genotype.toString());
			}
		}

		return result.toString();

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
		Set<GenotypeIDWithName> nonZeroGenotypes = new HashSet<GenotypeIDWithName>();
		for (GenotypeIDWithName genotypeIDWithName : genotypeIDsWithName)
		{
			Set<String> fragments = null;
			// MNV needs a bit of special treatment since only constituent SNVs
			// are stored in the pileup
			if (genotypeIDWithName.genotypeID.type == GenotypeEventType.MNV)
			{
				fragments = getMNVSupportingFragments(
						genotypeIDWithName.genotypeID);
			}
			else
			{
				fragments = pileup.genotypes.get(genotypeIDWithName.genotypeID);
			}

			if (fragments == null)
			{
				fragments = new HashSet<String>();
			}
			else
			{
				nonZeroGenotypes.add(genotypeIDWithName);
			}

			supportingFragments.put(genotypeIDWithName, fragments);
		}

		// compute powerset only for the genotypes with non-zero support
		Set<Set<GenotypeIDWithName>> powerSet = null;
		if (nonZeroGenotypes.size() <= 5)
		{
			// choose all non-zero genotypes
			powerSet = Sets.powerSet(nonZeroGenotypes);
		}
		else
		{
			// choose top 5 non-zero genotypes by number of supporting fragments
			Set<Entry<GenotypeIDWithName, Set<String>>> entries = supportingFragments
					.entrySet();
			List<Entry<GenotypeIDWithName, Set<String>>> list = new ArrayList<Entry<GenotypeIDWithName, Set<String>>>(
					entries);
			Collections.sort(list, fragmentCountComparator);

			Set<GenotypeIDWithName> chosenGenotypeIDs = new HashSet<GenotypeIDWithName>();
			for (int i = 0; i < 5; i++)
			{
				if (list.get(i).getValue().size() == 0)
				{
					break;
				}

				chosenGenotypeIDs.add(list.get(i).getKey());
			}

			powerSet = Sets.powerSet(chosenGenotypeIDs);
		}

		Set<Set<GenotypeIDWithName>> processingSet = new HashSet<Set<GenotypeIDWithName>>(
				powerSet);
		// add individual genotypes to the processing set
		for (GenotypeIDWithName genotypeIDWithName : genotypeIDsWithName)
		{
			Set<GenotypeIDWithName> s = new HashSet<GenotypeIDWithName>();
			s.add(genotypeIDWithName);
			processingSet.add(s);
		}

		powerSet = null;

		// populate the returning set
		Map<Set<GenotypeIDWithName>, Set<String>> returningSet = new HashMap<Set<GenotypeIDWithName>, Set<String>>();
		for (Set<GenotypeIDWithName> s : processingSet)
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

			// if any member genotype has 0 supporting fragments, return empty
			// set
			if (supporting.isEmpty())
			{
				return new HashSet<String>();
			}

			// first member
			// make sure to copy the set to avoid obscure bugs
			if (intersecting == null)
			{
				intersecting = new HashSet<String>(supporting);
			}
			else
			{
				intersecting.retainAll(supporting);
			}
		}

		return intersecting;
	}

	private Set<String> getMNVSupportingFragments(GenotypeID genotypeID)
	{
		List<GenotypeID> SNVs = MNVToSNVs(genotypeID);
		Set<String> intersecting = null;
		for (GenotypeID SNV : SNVs)
		{
			Set<String> fragments = pileup.genotypes.get(SNV);

			// if one of the constituent SNVs has zero support, then the MNV has
			// zero support
			if (fragments == null)
			{
				return null;
			}

			// first member
			// make sure to copy the set to avoid obscure bugs
			if (intersecting == null)
			{
				intersecting = new HashSet<String>(fragments);
			}
			else
			{
				intersecting.retainAll(fragments);
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
				continue;
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
