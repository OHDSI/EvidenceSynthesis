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

import dr.inference.operators.NormalGammaPrecisionGibbsOperator;

import dr.inference.distribution.GammaStatisticsProvider;
import org.ohdsi.metaAnalysis.ScalePrior.Base;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
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
//		return new NormalGammaPrecisionGibbsOperator(new GammaGibbsProvider.Default(hierarchicalLikelihood),
//				new NormalGammaPrecisionGibbsOperator.GammaParametrization(tauPrior.getDistribution()), null, weight);
		// updated line to use public re-implementation of `GammaParametrization`
		final int[] indices = {0};
		return new NormalGammaPrecisionGibbsOperator(new GammaGibbsProvider.Default(hierarchicalLikelihood),
				new GammaOnPrecisionProvider(tauPrior.getDistribution()), indices, weight);
//		return new NormalGammaPrecisionGibbsOperator(new GammaGibbsProvider.Default(hierarchicalLikelihood),
//				(GammaStatisticsProvider) tauPrior.getDistribution(), null, weight);
	}
}
