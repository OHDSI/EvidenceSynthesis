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

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

import dr.inference.distribution.EmpiricalDistributionData;
import dr.inference.distribution.EmpiricalDistributionLikelihood;
import dr.math.distributions.NormalDistribution;

public class SkewNormalDataModel extends DataModel.Base implements DataModel {

	public SkewNormalDataModel(String fileName) {
		this();
		File file = new File(fileName);
		List<List<String>> lines = new ArrayList<>();

		try {

			Scanner inputStream = new Scanner(file);

			while (inputStream.hasNext()) {
				String line = inputStream.next();
				line = line.replaceAll("\"", "");
				String[] values = line.split(",");
				lines.add(Arrays.asList(values));
			}

			inputStream.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		int muIdx = lines.get(0).indexOf("mu");
		int sigmaIdx = lines.get(0).indexOf("sigma");
		int gammaIdx = lines.get(0).indexOf("gamma");
		for (int i = 1; i < lines.size(); ++i) {
			double[] parameters = new double[] { Double.parseDouble(lines.get(i).get(muIdx)),
					Double.parseDouble(lines.get(i).get(sigmaIdx)), Double.parseDouble(lines.get(i).get(gammaIdx)) };
			addLikelihoodParameters(parameters, null);
		}
		finish();
	}

	public SkewNormalDataModel() {
	}

	@Override
	EmpiricalDistributionLikelihood makeFunctionalForm(List<EmpiricalDistributionData> dataList) {

		return new EmpiricalDistributionLikelihood(dataList, false) {

			private static final long serialVersionUID = 6515855145410583409L;
			
			private final double skewNormal(double x, double mean, double stdDev, double skew) {
				return Math.log(2) + NormalDistribution.logPdf(x, mean, stdDev) +
						NormalDistribution.cdf(skew * (x - mean), 0.0, stdDev, true);
			}

			@Override
			protected double logPDF(double x, EmpiricalDistributionData data) {

				final double mean = data.values[0];
				final double stdDev = data.values[1];
				final double skew = data.values[2];

				return skewNormal(x, mean, stdDev, skew);
			}

			@Override
			protected double gradientLogPdf(double x, EmpiricalDistributionData data) {

				final double mean = data.values[0];
				final double stdDev = data.values[1];
				final double skew = data.values[2];
				final double y = skew * (x - mean);

				// 2.  What about:
				return NormalDistribution.gradLogPdf(x, mean, stdDev) +
						skew * NormalDistribution.pdf(y, 0.0, stdDev) /
								NormalDistribution.cdf(y, 0.0, stdDev, false);
			}
		};
	}
}
