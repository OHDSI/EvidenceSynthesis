package org.ohdsi.metaAnalysis;

import org.ohdsi.metaAnalysis.ScalePrior.Base;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.NormalGammaPrecisionGibbsOperator;
import dr.inference.operators.repeatedMeasures.GammaGibbsProvider;
import dr.math.distributions.GammaDistribution;

public class GammaOnPrecisionPrior extends Base implements ScalePrior {

	public GammaOnPrecisionPrior(double scale, double rate) {
		super(new GammaDistribution(scale, rate));
	}

	@Override
	public boolean isPrecision() {
		return true;
	}

	@Override
	public MCMCOperator getOperator(DistributionLikelihood hierarchicalLikelihood, double weight, AdaptationMode mode) {
		return new NormalGammaPrecisionGibbsOperator(new GammaGibbsProvider.Default(hierarchicalLikelihood),
				tauPrior.getDistribution(), null, weight);
	}
}
