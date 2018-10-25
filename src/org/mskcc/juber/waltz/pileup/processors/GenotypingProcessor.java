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
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
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

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;

public class GenotypingProcessor implements PileupProcessor
{
	private RegionPileupView pileup;
	private PileupMetricsProcessor metricsProcessor;
	private Set<GenotypeIDWithMafLine> genotypeIDsWithMafLine;
	private Comparator<Entry<GenotypeIDWithMafLine, Set<String>>> fragmentCountComparator;
	private String mafHeader;
	private Map<String, Integer> mafColumns;
	private IndexedFastaSequenceFile referenceFasta;

	public GenotypingProcessor(File lociFile,
			IndexedFastaSequenceFile referenceFasta) throws IOException
	{
		this.referenceFasta = referenceFasta;

		genotypeIDsWithMafLine = new HashSet<GenotypeIDWithMafLine>();

		BufferedReader reader = new BufferedReader(new FileReader(lociFile));
		String line = null;
		String[] parts = null;

		processMafHeader(reader.readLine());

		while ((line = reader.readLine()) != null)
		{
			parts = line.split("\t");

			// expanded parts include additional Waltz fields
			String[] expandedParts = new String[mafColumns.size()];

			// copy the original parts into expanded parts
			for (int i = 0; i < expandedParts.length; i++)
			{
				if (i < parts.length)
				{
					expandedParts[i] = parts[i];
				}

				// remove nulls
				if (expandedParts[i] == null)
				{
					expandedParts[i] = "";
				}
			}

			GenotypeID genotypeID = makeGenotypeID(expandedParts);
			if (genotypeID == null)
			{
				continue;
			}

			genotypeIDsWithMafLine
					.add(new GenotypeIDWithMafLine(genotypeID, expandedParts));
		}

		reader.close();

		// create fragment count comparator
		fragmentCountComparator = new Comparator<Entry<GenotypeIDWithMafLine, Set<String>>>()
		{

			@Override
			public int compare(Entry<GenotypeIDWithMafLine, Set<String>> arg0,
					Entry<GenotypeIDWithMafLine, Set<String>> arg1)
			{
				return arg1.getValue().size() - arg0.getValue().size();
			}
		};
	}

	private void processMafHeader(String header)
	{
		mafColumns = new LinkedHashMap<String, Integer>();

		// Add Waltz fields
		mafHeader = header + "\tWaltz_total_t_depth"
				+ "\tWaltz_total_t_alt_count" + "\tWaltz_MD_t_depth"
				+ "\tWaltz_MD_t_alt_count";

		String[] parts = mafHeader.split("\t");

		for (int i = 0; i < parts.length; i++)
		{
			mafColumns.put(parts[i], i);
		}
	}

	public String getMafHeader()
	{
		return mafHeader;
	}

	public Set<Interval> getGenotypesAsIntervals()
	{
		Set<Interval> intervals = new HashSet<Interval>();

		for (GenotypeIDWithMafLine genotypeIDWithMafLine : genotypeIDsWithMafLine)
		{
			if (genotypeIDWithMafLine == null
					|| genotypeIDWithMafLine.genotypeID == null)
			{
				continue;
			}

			GenotypeID id = genotypeIDWithMafLine.genotypeID;
			intervals.add(new Interval(id.contig, id.position, id.endPosition,
					false, genotypeIDWithMafLine.name));
		}

		return intervals;
	}

	/**
	 * make genotype ids from a parsed line from a maf file
	 * 
	 * @param parts
	 * @return
	 */
	private GenotypeID makeGenotypeID(String[] parts)
	{
		String contig = parts[mafColumns.get("Chromosome")];
		int position = Integer
				.parseInt(parts[mafColumns.get("Start_Position")]);
		String refString = parts[mafColumns.get("Reference_Allele")];
		String altString = parts[mafColumns.get("Tumor_Seq_Allele2")];
		String eventTypeString = parts[mafColumns.get("Variant_Type")];

		// complex event
		if ((eventTypeString.equals("INS") || eventTypeString.equals("DEL"))
				&& !refString.equals("-") && !altString.equals("-"))
		{
			String line = StringUtil.join("\t", Arrays.asList(parts));

			System.out.println("Complex events not supported yet!");
			System.out.println(line);
			return null;
		}

		if (eventTypeString.equals("SNP"))
		{
			byte[] ref = new byte[] { (byte) refString.charAt(0) };
			byte[] alt = new byte[] { (byte) altString.charAt(0) };

			return new GenotypeID(GenotypeEventType.SNV, contig, position, ref,
					alt);
		}
		else if (eventTypeString.equals("DNP") || eventTypeString.equals("TNP")
				|| eventTypeString.equals("MNP")
				|| eventTypeString.equals("MNV")
				|| eventTypeString.equals("ONP"))
		{
			byte[] ref = new byte[refString.length()];
			for (int i = 0; i < ref.length; i++)
			{
				ref[i] = (byte) refString.charAt(i);
			}

			byte[] alt = new byte[altString.length()];
			for (int i = 0; i < alt.length; i++)
			{
				alt[i] = (byte) altString.charAt(i);
			}

			return new GenotypeID(GenotypeEventType.MNV, contig, position, ref,
					alt);
		}
		else if (eventTypeString.equals("INS"))
		{
			byte[] ref = referenceFasta
					.getSubsequenceAt(contig, position, position).getBases();
			byte[] alt = new byte[altString.length() + 1];
			alt[0] = ref[0];
			for (int i = 1; i < alt.length; i++)
			{
				alt[i] = (byte) altString.charAt(i - 1);
			}

			return new GenotypeID(GenotypeEventType.INSERTION, contig, position,
					ref, alt);
		}
		else if (eventTypeString.equals("DEL"))
		{
			position--;
			byte[] alt = referenceFasta
					.getSubsequenceAt(contig, position, position).getBases();
			byte[] ref = new byte[refString.length() + 1];
			ref[0] = alt[0];
			for (int i = 1; i < ref.length; i++)
			{
				ref[i] = (byte) refString.charAt(i - 1);
			}

			return new GenotypeID(GenotypeEventType.DELETION, contig, position,
					ref, alt);
		}
		else
		{
			System.err.println("Not Supported: " + contig + "\t" + position
					+ "\t" + eventTypeString + "\t" + refString + "\t"
					+ altString + "\n");
			return null;
		}
	}

	/**
	 * make genotype ids from a parsed line from an input mutations files
	 * 
	 * @param parts
	 * @return
	 */
	private GenotypeID makeGenotypeIDOld(String[] parts)
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
		Set<GenotypeIDWithMafLine> regionGenotypes = new HashSet<GenotypeIDWithMafLine>();

		// pick out the genotypes contained in the current region/pileup
		// for (int i = 0; i < genotypeIDsWithMafLine.size(); i++)
		for (GenotypeIDWithMafLine genotypeIDWithMafLine : genotypeIDsWithMafLine)
		{
			if (genotypeIDWithMafLine == null
					|| genotypeIDWithMafLine.genotypeID == null)
			{
				continue;
			}

			if (pileup.contains(genotypeIDWithMafLine.genotypeID))
			{
				regionGenotypes.add(genotypeIDWithMafLine);
			}
		}

		// no genotypes to be profiled in the current region
		if (regionGenotypes.isEmpty())
		{
			return;
		}

		String outString = processGenotypeSet(regionGenotypes);
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
	private String processGenotypeSet(
			Set<GenotypeIDWithMafLine> regionGenotypes)
	{
		StringBuilder result = null;

		Map<Set<GenotypeIDWithMafLine>, Set<String>> powerSetFragments = getSupportingFragmentsForPowerSet(
				regionGenotypes);

		// process each element of the powerset
		for (Set<GenotypeIDWithMafLine> s : powerSetFragments.keySet())
		{
			Genotype genotype = null;

			Iterator<GenotypeIDWithMafLine> it = s.iterator();
			GenotypeIDWithMafLine firstGenotypeIDWithMafLine = it.next();
			String[] mafLineParts = firstGenotypeIDWithMafLine.mafLineParts
					.clone();

			// a single, traditional genotype
			if (s.size() == 1)
			{
				genotype = new Genotype(firstGenotypeIDWithMafLine.name,
						firstGenotypeIDWithMafLine.genotypeID);
			}
			else
			{
				StringBuilder name = null;
				// build name
				for (GenotypeIDWithMafLine genotypeIDWitMafLine : s)
				{
					if (name == null)
					{
						name = new StringBuilder(genotypeIDWitMafLine.name);
					}
					else
					{
						name.append("+");
						name.append(genotypeIDWitMafLine.name);
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
				result = new StringBuilder(
						makeOutputMafLine(genotype, mafLineParts));
			}
			else
			{
				result.append(System.lineSeparator());
				result.append(makeOutputMafLine(genotype, mafLineParts));
			}
		}

		return result.toString();

	}

	/**
	 * make the line that will go to the output as a maf line
	 * 
	 * @param genotype
	 * @param mafLineParts
	 * @return
	 */
	private String makeOutputMafLine(Genotype genotype, String[] mafLineParts)
	{
		// composite genotype, delete all fields and add composite name to
		// Variant_Type
		if (genotype.id == null)
		{
			// clear all fields
			for (int i = 0; i < mafLineParts.length; i++)
			{
				mafLineParts[i] = "";
			}

			// add composite name
			mafLineParts[mafColumns.get("Variant_Type")] = "COMPOSITE:"
					+ genotype.name;
		}

		// add Waltz genotyping info
		mafLineParts[mafColumns.get("Waltz_total_t_depth")] = ""
				+ genotype.totalCoverage;
		mafLineParts[mafColumns.get("Waltz_total_t_alt_count")] = ""
				+ genotype.totalSupportingCoverage;
		mafLineParts[mafColumns.get("Waltz_MD_t_depth")] = ""
				+ genotype.uniqueCoverage;
		mafLineParts[mafColumns.get("Waltz_MD_t_alt_count")] = ""
				+ genotype.uniqueSupportingCoverage;

		// make a string
		StringBuilder s = new StringBuilder(mafLineParts[0]);
		for (int i = 1; i < mafLineParts.length; i++)
		{
			s.append("\t");
			s.append(mafLineParts[i]);
		}

		return s.toString();
	}

	/**
	 * Iterate over the powerset of the given set of genotypes and find the
	 * fragments supporting each non-empty element of the powerset
	 * 
	 * @param sNVs
	 * @return
	 */
	private Map<Set<GenotypeIDWithMafLine>, Set<String>> getSupportingFragmentsForPowerSet(
			Set<GenotypeIDWithMafLine> genotypeIDsWithName)
	{
		// find supporting fragments for individual genotypes
		Map<GenotypeIDWithMafLine, Set<String>> supportingFragments = new HashMap<GenotypeIDWithMafLine, Set<String>>();
		Set<GenotypeIDWithMafLine> nonZeroGenotypes = new HashSet<GenotypeIDWithMafLine>();
		for (GenotypeIDWithMafLine genotypeIDWithName : genotypeIDsWithName)
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
		Set<Set<GenotypeIDWithMafLine>> powerSet = null;
		if (nonZeroGenotypes.size() <= 5)
		{
			// choose all non-zero genotypes
			powerSet = Sets.powerSet(nonZeroGenotypes);
		}
		else
		{
			// choose top 5 non-zero genotypes by number of supporting fragments
			Set<Entry<GenotypeIDWithMafLine, Set<String>>> entries = supportingFragments
					.entrySet();
			List<Entry<GenotypeIDWithMafLine, Set<String>>> list = new ArrayList<Entry<GenotypeIDWithMafLine, Set<String>>>(
					entries);
			Collections.sort(list, fragmentCountComparator);

			Set<GenotypeIDWithMafLine> chosenGenotypeIDs = new HashSet<GenotypeIDWithMafLine>();
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

		Set<Set<GenotypeIDWithMafLine>> processingSet = new HashSet<Set<GenotypeIDWithMafLine>>(
				powerSet);
		// add individual genotypes to the processing set
		for (GenotypeIDWithMafLine genotypeIDWithName : genotypeIDsWithName)
		{
			Set<GenotypeIDWithMafLine> s = new HashSet<GenotypeIDWithMafLine>();
			s.add(genotypeIDWithName);
			processingSet.add(s);
		}

		powerSet = null;

		// populate the returning set
		Map<Set<GenotypeIDWithMafLine>, Set<String>> returningSet = new HashMap<Set<GenotypeIDWithMafLine>, Set<String>>();
		for (Set<GenotypeIDWithMafLine> s : processingSet)
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
	private Set<String> getIntersection(Set<GenotypeIDWithMafLine> genotypeIDs,
			Map<GenotypeIDWithMafLine, Set<String>> supportingFragments)
	{
		Set<String> intersecting = null;
		for (GenotypeIDWithMafLine genotypeIDWithName : genotypeIDs)
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
	private ContineousSpan computeSpan(Set<GenotypeIDWithMafLine> s)
	{
		String contig = null;
		int start = 0;
		int end = 0;

		for (GenotypeIDWithMafLine genotypeIDWithName : s)
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

	private class GenotypeIDWithMafLine
	{
		public final GenotypeID genotypeID;
		public final String name;
		String[] mafLineParts;

		public GenotypeIDWithMafLine(GenotypeID genotypeID,
				String[] mafLineParts)
		{
			this.genotypeID = genotypeID;
			this.mafLineParts = mafLineParts;
			this.name = makeGenotypeName(mafLineParts);

			clearSampleSpecificFields();
		}

		private String makeGenotypeName(String[] mafLineParts)
		{
			String name = null;
			if (mafColumns.get("Hugo_Symbol") == null
					|| mafColumns.get("HGVSp_Short") == null
					|| mafLineParts[mafColumns.get("Hugo_Symbol")].equals("")
					|| mafLineParts[mafColumns.get("HGVSp_Short")].equals(""))
			{
				String contig = mafLineParts[mafColumns.get("Chromosome")];
				int position = Integer.parseInt(
						mafLineParts[mafColumns.get("Start_Position")]);
				String refString = mafLineParts[mafColumns
						.get("Reference_Allele")];
				String altString = mafLineParts[mafColumns
						.get("Tumor_Seq_Allele2")];
				String eventTypeString = mafLineParts[mafColumns
						.get("Variant_Type")];

				name = contig + "-" + position + "-" + eventTypeString + "-"
						+ refString + "-" + altString;
			}
			else
			{
				String proteinChange = mafLineParts[mafColumns
						.get("HGVSp_Short")];
				if (proteinChange.startsWith("p."))
				{
					proteinChange = proteinChange.substring(2);
				}

				name = mafLineParts[mafColumns.get("Hugo_Symbol")] + "."
						+ proteinChange;
			}

			return name;

		}

		/**
		 * from the maf line parts, clear the sample-specific fields that don't
		 * make sense in the context of genotyping
		 */
		private void clearSampleSpecificFields()
		{
			String[] sampleSpecificFields = new String[] {
					"Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
					"Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
					"Tumor_Validation_Allele1", "Tumor_Validation_Allele2",
					"Match_Norm_Validation_Allele1",
					"Match_Norm_Validation_Allele2", "Verification_Status",
					"Validation_Status", "Mutation_Status", "Sequencing_Phase",
					"Sequence_Source", "Validation_Method", "Score", "BAM_File",
					"Sequencer", "Tumor_Sample_UUID",
					"Matched_Norm_Sample_UUID", "t_depth", "t_ref_count",
					"t_alt_count", "n_depth", "n_ref_count", "n_alt_count" };

			for (int i = 0; i < sampleSpecificFields.length; i++)
			{
				Integer columnNumber = mafColumns.get(sampleSpecificFields[i]);
				if (columnNumber == null)
				{
					continue;
				}

				mafLineParts[columnNumber] = "";
			}

		}

		@Override
		public int hashCode()
		{
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((genotypeID == null) ? 0 : genotypeID.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj)
		{
			if (this == obj)
			{
				return true;
			}
			if (obj == null)
			{
				return false;
			}
			if (getClass() != obj.getClass())
			{
				return false;
			}
			GenotypeIDWithMafLine other = (GenotypeIDWithMafLine) obj;

			if (genotypeID == null)
			{
				if (other.genotypeID != null)
				{
					return false;
				}
			}
			else if (!genotypeID.equals(other.genotypeID))
			{
				return false;
			}

			return true;
		}
	}
}
