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
package org.mskcc.juber.waltz.pileup.processors.signatures;

/**
 * @author Juber Patel
 * 
 *         A locus that matches a signature.
 * 
 *         A simple bean class
 *
 */
public class SignatureLocus
{
	public final int start;
	public final int end;
	public final String description;
	public final String evidence;

	public SignatureLocus(int start, int end, String description,
			String evidence)
	{
		this.start = start;
		this.end = end;
		this.description = description;
		this.evidence = evidence;
	}
}
