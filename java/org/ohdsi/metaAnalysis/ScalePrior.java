package org.ohdsi.metaAnalysis;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.distributions.Distribution;

public interface ScalePrior {

	Parameter getParameter();

	boolean isPrecision();

	Likelihood getPrior();

	MCMCOperator getOperator(DistributionLikelihood hierarchicalLikelihood, double weight, AdaptationMode mode);

	abstract class Base implements ScalePrior {

		final Parameter tau;
		final DistributionLikelihood tauPrior;

		protected Base(Distribution distribution) {
			tau = new Parameter.Default("tau", 2.0, 0.0, Double.POSITIVE_INFINITY);
			tauPrior = new DistributionLikelihood(distribution);
			tauPrior.addData(tau);
		}

		@Override
		public Parameter getParameter() {
			return tau;
		}

		@Override
		public Likelihood getPrior() {
			return tauPrior;
		}
	}
}
