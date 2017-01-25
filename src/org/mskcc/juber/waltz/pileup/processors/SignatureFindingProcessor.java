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
package org.mskcc.juber.waltz.pileup.processors;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.mskcc.juber.waltz.WaltzOutput;
import org.mskcc.juber.waltz.pileup.RegionPileupView;
import org.mskcc.juber.waltz.pileup.processors.signatures.PileupSignature;
import org.mskcc.juber.waltz.pileup.processors.signatures.SignatureLocus;
import org.mskcc.juber.waltz.pileup.processors.signatures.TranslocationBreakpointSignature;

/**
 * @author Juber Patel
 *
 */
public class SignatureFindingProcessor implements PileupProcessor
{
	private RegionPileupView pileup;
	private List<PileupSignature> signatures;

	public SignatureFindingProcessor(String moduleArgument)
	{
		setSignaturesToLookFor(moduleArgument);
	}

	@Override
	public void setRegionPileupView(RegionPileupView view)
	{
		this.pileup = view;
	}

	@Override
	public void processRegion(WaltzOutput output) throws IOException
	{
		String sampleName = output.getSampleName();
		String contig = pileup.interval.getContig();

		for (int i = 0; i < signatures.size(); i++)
		{
			PileupSignature signature = signatures.get(i);
			List<SignatureLocus> loci = signature.findIn(pileup);

			// for each locus that fits the current signature
			for (SignatureLocus locus : loci)
			{
				String outString = sampleName + "\t" + contig + "\t"
						+ locus.start + "\t" + locus.end + "\t"
						+ locus.description + "\t" + locus.evidence + "\n";
				output.toSignatureIntervalsWriter(outString);
			}
		}

		// free the memory once we are done
		this.pileup = null;
	}

	private void setSignaturesToLookFor(String moduleArgument)
	{
		signatures = new ArrayList<PileupSignature>();

		String[] signaturesRequested = moduleArgument.split(",");
		for (String sig : signaturesRequested)
		{
			if (sig.equals("TranslocationBreakpoint"))
			{
				signatures.add(new TranslocationBreakpointSignature());
			}
		}

	}

}
