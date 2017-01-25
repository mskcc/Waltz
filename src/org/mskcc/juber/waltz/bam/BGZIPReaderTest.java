/*******************************************************************************
 * @author Juber Patel
 * 
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
package org.mskcc.juber.waltz.bam;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

import htsjdk.samtools.util.BlockCompressedInputStream;

/**
 * @author Juber Patel
 * 
 */
public class BGZIPReaderTest
{

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{
		// this is actually slower than SamReader. 15 seconds vs 12 seconds.
		// also, not printing all characters properly.
		// this may be because BinaryCodec forces little-endian reading and
		// writing.
		// GZipInputStream is even slower: 19 seconds

		BufferedReader reader = new BufferedReader(
				new InputStreamReader(new BlockCompressedInputStream(new File(
						"SSC12236-1-small.bam"))));

		String line = null;
		int count = 0;
		long start = System.currentTimeMillis();

		while ((line = reader.readLine()) != null)
		{
			// count++;
			// if(count == 300)
			// break;
			//
			// System.out.println(line);

		}

		long time = System.currentTimeMillis() - start;

		System.out.println((time * 1.0) / 1000);

	}
}
