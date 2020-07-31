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

				if (x <= data.values[0]) {
					return data.density[0];
				} else if (x >= data.values[end]) {
					return data.density[end];
				} else {
					return super.logPDF(x, data);
				}
			}

			@Override
			public double gradientLogPdf(double x, EmpiricalDistributionData data) {
				final int end = data.values.length;

				if (x <= data.values[0] || x >= data.values[end]) {
					return 0.0;
				} else {
					return super.gradientLogPdf(x, data);
				}
			}
		};
	}
}
