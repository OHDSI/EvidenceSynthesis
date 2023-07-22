/*******************************************************************************
 * Copyright 2023 Observational Health Data Sciences and Informatics
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

import dr.inference.distribution.EmpiricalDistributionData;
import dr.inference.distribution.EmpiricalDistributionLikelihood;
import dr.inference.distribution.SplineInterpolatedLikelihood;

import java.util.List;

public class ExtendingEmpiricalDataModel extends EmpiricalDataModel {

	public ExtendingEmpiricalDataModel(String fileName) {
		super(fileName);
	}

	public ExtendingEmpiricalDataModel() {
	}

	@Override
	EmpiricalDistributionLikelihood makeFunctionalForm(List<EmpiricalDistributionData> dataList) {
		return new SplineInterpolatedLikelihood(dataList, 1, false) {

			private static final long serialVersionUID = 6586114515564784388L;

			@Override
			public double logPDF(double x, EmpiricalDistributionData data) {
				final int end = data.values.length - 1;
				// Assume density is constant outside of range:
				//				if (x <= data.values[0]) {
				//					return data.density[0];
				//				} else if (x >= data.values[end]) {
				//					return data.density[end];
				//				} else {
				//					return super.logPDF(x, data);
				//				}

				// Linear extrapolation outside of range:
				if (end == 0) {
					return data.density[0];
				} else if (x <= data.values[0]) {
					final double slope = Math.max(0, (data.density[1] - data.density[0]) / (data.values[1] - data.values[0]));
					return (x - data.values[0]) * slope + data.density[0];
				} else if (x >= data.values[end]) {
					final double slope = Math.min(0, (data.density[end] - data.density[end - 1]) / (data.values[end] - data.values[end - 1]));
					return (x - data.values[end]) * slope + data.density[end];
				} else {
					return super.logPDF(x, data);
				}

				// Second-order extrapolation outside of range:
				//				if (x <= data.values[0]) {
				//					final double slope1 = (data.density[1] - data.density[0]) / (data.values[1] - data.values[0]);
				//					final double slope2 = (data.density[2] - data.density[1]) / (data.values[2] - data.values[1]);
				//					final double secondD = (slope2 - slope1) / (((data.values[2] + data.values[1]) / 2.0) - ((data.values[1] + data.values[0]) / 2.0));
				//					return data.density[0] + slope1 * (x - data.values[0]) + 0.5 * secondD * sqr(x - data.values[0]);
				//				} else if (x >= data.values[end]) {
				//					final double slope1 = (data.density[end-1] - data.density[end-2]) / (data.values[end-1] - data.values[end-2]);
				//					final double slope2 = (data.density[end] - data.density[end-1]) / (data.values[end] - data.values[end-1]);
				//					final double secondD = (slope2 - slope1) / (((data.values[end] + data.values[end-1]) / 2.0) - ((data.values[end-1] + data.values[end-2]) / 2.0));
				//					return data.density[end] + slope2 * (x - data.values[end]) + 0.5 * secondD * sqr(x - data.values[end]);
				//				} else {
				//					return super.logPDF(x, data);
				//				}
			}

			public final double sqr(double x) {
				return x * x;
			}

			@Override
			public double gradientLogPdf(double x, EmpiricalDistributionData data) {
				final int end = data.values.length;

				if (x <= data.values[0] || x >= data.values[end - 1]) {
					return 0.0;
				} else {
					return super.gradientLogPdf(x, data);
				}
			}
		};
	}
}
