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
