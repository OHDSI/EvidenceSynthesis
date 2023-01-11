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
