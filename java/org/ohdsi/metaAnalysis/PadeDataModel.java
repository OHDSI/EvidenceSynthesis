/*******************************************************************************
 * Copyright 2021 Observational Health Data Sciences and Informatics
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
package org.ohdsi.metaAnalysis;

import java.util.List;

import dr.inference.distribution.EmpiricalDistributionData;
import dr.inference.distribution.EmpiricalDistributionLikelihood;

public class PadeDataModel extends DataModel.Base implements DataModel {
	
	public PadeDataModel() {
	}
	
	@Override
	EmpiricalDistributionLikelihood makeFunctionalForm(List<EmpiricalDistributionData> dataList) {
		
		return new EmpiricalDistributionLikelihood(dataList, false) {
			
			private static final long serialVersionUID = 6515855145410583409L;
			
			private final double sqr(double x) {
				return x*x;
			}
			
			private final double pade(double x, double beta, double a0, double a1, double a2, double b1, double b2) {
				final double delta = x - beta;
				final double numerator = a0 + delta * a1 + sqr(delta) * a2;
				final double denominator = 1 + delta * b1 + sqr(delta) * b2;
				return numerator / denominator;
			}
			
			private final double padeGradient(double x, double beta, double a0, double a1, double a2, double b1, double b2) {
				final double delta = x - beta;
				final double p = a0 + delta * a1 + sqr(delta) * a2;
				final double q = 1 + delta * b1 + sqr(delta) * b2;
				final double gradientNumerator =  q*(a1 + 2*a2*delta) - p*(b1+2*b2*delta);
				return gradientNumerator / (sqr(q));
			}
			
			@Override
			protected double logPDF(double x, EmpiricalDistributionData data) {
				final double beta = data.values[0];
				final double a0 = data.values[1];
				final double a1 = data.values[2];
				final double a2 = data.values[3];
				final double b1 = data.values[4];
				final double b2 = data.values[5];
				final double minBeta = data.values[6];
				final double maxBeta = data.values[7];
				final double minD1 = data.values[8];
				final double maxD1 = data.values[9];
				if (x < minBeta) 
					return pade(minBeta, beta, a0, a1, a2, b1, b2) + (x - minBeta) * minD1;
				else if (x > maxBeta)
					return pade(maxBeta, beta, a0, a1, a2, b1, b2) + (x - maxBeta) * maxD1;
				else	
					return pade(x, beta, a0, a1, a2, b1, b2);
			}
			
			@Override
			protected double gradientLogPdf(double x, EmpiricalDistributionData data) {
				final double beta = data.values[0];
				final double a0 = data.values[1];
				final double a1 = data.values[2];
				final double a2 = data.values[3];
				final double b1 = data.values[4];
				final double b2 = data.values[5];
				final double minBeta = data.values[6];
				final double maxBeta = data.values[7];
				final double minD1 = data.values[8];
				final double maxD1 = data.values[9];
				if (x < minBeta) 
					return minD1;
				else if (x > maxBeta)
					return maxD1;
				else	
					return padeGradient(x, beta, a0, a1, a2, b1, b2);
			}
		};
	}
}
