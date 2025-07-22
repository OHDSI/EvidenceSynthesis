/*******************************************************************************
 * Copyright 2025 Observational Health Data Sciences and Informatics
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
import dr.inference.distribution.SplineInterpolatedLikelihood;

public class EmpiricalDataModel extends DataModel.Base implements DataModel {

	public EmpiricalDataModel(String fileName) {
		this();
		File file = new File(fileName);
		List<List<String>> lines = new ArrayList<>();

		try {

			Scanner inputStream = new Scanner(file);

			while (inputStream.hasNext()) {
				String line = inputStream.next();
				String[] values = line.split(",");
				lines.add(Arrays.asList(values));
			}

			inputStream.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		double[] x = parseRow(lines.get(0), true);

		for (int i = 1; i < lines.size(); ++i)
			addLikelihoodParameters(x, parseRow(lines.get(i), false));

		finish();
	}

	public EmpiricalDataModel() {
	}

	private static double[] parseRow(List<String> row, boolean strip) {
		double[] values = new double[row.size()];
		for (int i = 0; i < row.size(); ++i) {
			String string = row.get(i);
			if (strip) {
				string = string.replace("\"", "");
			}
			values[i] = Double.parseDouble(string);
		}
		return values;
	}

	@Override
	EmpiricalDistributionLikelihood makeFunctionalForm(List<EmpiricalDistributionData> dataList) {
		return new SplineInterpolatedLikelihood(dataList, 1, false);
	}
}
