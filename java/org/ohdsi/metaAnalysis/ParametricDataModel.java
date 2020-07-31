package org.ohdsi.metaAnalysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

import dr.inference.distribution.EmpiricalDistributionData;
import dr.inference.distribution.EmpiricalDistributionLikelihood;

public class ParametricDataModel extends DataModel.Base implements DataModel {

	public ParametricDataModel(String fileName) {
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

	public ParametricDataModel() {
	}

	@Override
	EmpiricalDistributionLikelihood makeFunctionalForm(List<EmpiricalDistributionData> dataList) {

		return new EmpiricalDistributionLikelihood(dataList, false) {

			private static final long serialVersionUID = 5902404628932924630L;

			private final double sqr(double x) {
				return x * x;
			}

			@Override
			protected double logPDF(double x, EmpiricalDistributionData data) {
				final double mu = data.values[0];
				final double sigma = data.values[1];
				final double gamma = data.values[2];

				return ((Math.exp(gamma * (x - mu)))) * ((-sqr(x - mu)) / (2 * sqr(sigma)));
			}

			@Override
			protected double gradientLogPdf(double x, EmpiricalDistributionData data) {
				final double mu = data.values[0];
				final double sigma = data.values[1];
				final double gamma = data.values[2];

				return -(Math.exp(gamma * (x - mu)) * (gamma * (x - mu) + 2) * (x - mu)) / (2 * sqr(sigma));
			}
		};
	}
}
