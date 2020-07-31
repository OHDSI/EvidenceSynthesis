package org.ohdsi.metaAnalysis;

import org.ohdsi.metaAnalysis.ScalePrior.Base;

import dr.inference.distribution.CauchyDistribution;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.ScaleOperator;

public class HalfCauchyOnStdDevPrior extends Base implements ScalePrior {

	public HalfCauchyOnStdDevPrior(double median, double scale) {
		super(new CauchyDistribution(median, scale));
	}

	@Override
	public boolean isPrecision() {
		return false;
	}

	@Override
	public MCMCOperator getOperator(DistributionLikelihood hierarchicalLikelihood, double weight, AdaptationMode mode) {
		return new ScaleOperator(tau, 0.75, mode, weight);
	}
}
