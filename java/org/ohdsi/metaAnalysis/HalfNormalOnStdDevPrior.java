package org.ohdsi.metaAnalysis;

import org.ohdsi.metaAnalysis.ScalePrior.Base;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.ScaleOperator;
import dr.math.distributions.NormalDistribution;

public class HalfNormalOnStdDevPrior extends Base implements ScalePrior {

	public HalfNormalOnStdDevPrior(double mean, double standardDeviation) {
		super(new NormalDistribution(mean, standardDeviation));
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
