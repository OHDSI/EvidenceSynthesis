/*******************************************************************************
 * Copyright 2024 Observational Health Data Sciences and Informatics
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package org.ohdsi.data;

/**
 * @author Martijn Schuemie
 */
public class SccsData {
	public int[] y;
	public double[] a;
	public double[][] x;
	public int[] stratumId;
	public double[] time;

	public SccsData(int[] y, double[] a, double[][] x, int[] stratumId, double[] time) {
		if (y.length < 2)
			throw new IllegalArgumentException("SCCS data must have at least two rows");
		if (y.length != time.length || y.length != x.length || y.length != a.length || y.length != stratumId.length) 
			throw new IllegalArgumentException("All dimensions must be equal");
		for (int i = 1; i < stratumId.length; i++) 
			if (stratumId[i-1] > stratumId[i]) 
				throw new IllegalArgumentException("SCCS data must be sorted by stratumId");
		this.y = y;
		this.a = a;
		this.x = x;
		this.stratumId = stratumId;
		this.time = time;
	}
}