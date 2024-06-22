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
package org.ohdsi.likelihood;

import dr.inference.model.*;
import org.ohdsi.data.SortedCoxData;

/**
 * @author Marc A. Suchard
 */
public class MultivariableCoxPartialLikelihood extends ConditionalPoissonLikelihood {

	public MultivariableCoxPartialLikelihood(Parameter beta, Parameter nuisance, SortedCoxData data) {
		super("multivariableCoxPartialLikelihood", beta, nuisance, data);
	}

	@Override
	protected double calculateLogLikelihood() {

		final int[] y = data.y;
		final double[] x = data.x;
		final double[] weight = data.weight;
		final int[] strata = data.strata;

		final double[] beta = this.beta.getParameterValues(); // TODO Faster as copy or direct access?

		double logLikelihood = 0.0;
		int resetIndex = 0;
		double denominator = 0.0;

		for (int i = 0; i < N; ++i) {
			if (i == strata[resetIndex]) {
				denominator = 0.0;
				++resetIndex;
			}

			double rowXBeta = 0.0;
			for (int j = 0; j < J; ++j) {
				rowXBeta += x[j * N + i] * beta[j]; // TODO Need to form xTranspose
			}

			denominator += weight[i] * Math.exp(rowXBeta);
			logLikelihood += y[i] * weight[i] * (rowXBeta - Math.log(denominator)); // TODO Could pre-compute total numerator as function of beta
		}

		return logLikelihood;
	}

	public static void main(String[] args) {

		int[] y;
		double[] x;
		double[] weight;
		int[] strata;
		double beta;

		SortedCoxData data;
		Parameter parameter;
		Likelihood cox;

		// No ties, no strata, no weights
		y = new int[] { 1, 1, 0, 1, 1, 0, 1 };
		x = new double[] { 0, 2, 0, 0, 1, 1, 1 };
		weight = new double[] { 1, 1, 1, 1, 1, 1, 1 };
		strata = new int[] { 7 };
		beta = 0.3883064;
		// logLik = -5.401371

		data = new SortedCoxData(y, x, strata, weight);
		parameter = new Parameter.Default(beta);
		cox = new MultivariableCoxPartialLikelihood(parameter, null, data);
		System.err.println("M " + cox.getLogLikelihood());
		data = new SortedCoxData(y, x, strata, weight);
		parameter = new Parameter.Default(beta);
		cox = new CoxPartialLikelihood(parameter, data);
		System.err.println("C " + cox.getLogLikelihood());

		System.err.println();

		// No ties, no strata, with weights
		weight = new double[] { 1, 1, 1, 1, 0, 1, 1 };
		beta = 0.3443102;
		// logLik = -3.726085

		data = new SortedCoxData(y, x, strata, weight);
		parameter = new Parameter.Default(beta);
		cox = new MultivariableCoxPartialLikelihood(parameter, null, data);
		System.err.println("M " + cox.getLogLikelihood());
		data = new SortedCoxData(y, x, strata, weight);
		parameter = new Parameter.Default(beta);
		cox = new CoxPartialLikelihood(parameter, data);
		System.err.println("C " + cox.getLogLikelihood());

		System.err.println();

		// No ties, with strata, no weights
		y = new int[] { 1, 1, 0, 1, 0, 1, 1 };
		weight = new double[] { 1, 1, 1, 1, 1, 1, 1 };
		x = new double[] { 0, 2, 1, 1, 0, 0, 1 };
		strata = new int[] { 4, 7 };
		beta = 1.205852;
		// logLike = -2.978028

		data = new SortedCoxData(y, x, strata, weight);
		parameter = new Parameter.Default(beta);
		cox = new MultivariableCoxPartialLikelihood(parameter, null, data);
		System.err.println("M " + cox.getLogLikelihood());
		data = new SortedCoxData(y, x, strata, weight);
		parameter = new Parameter.Default(beta);
		cox = new CoxPartialLikelihood(parameter, data);
		System.err.println("C " + cox.getLogLikelihood());

		System.err.println();

		// No ties, with strata, with weights
		weight = new double[] { 1, 1, 1, 1, 1, 1, 0 };
		beta = 0.7563076;
		// logLike = -2.418282

		data = new SortedCoxData(y, x, strata, weight);
		parameter = new Parameter.Default(beta);
		cox = new MultivariableCoxPartialLikelihood(parameter, null, data);
		System.err.println("M " + cox.getLogLikelihood());
		data = new SortedCoxData(y, x, strata, weight);
		parameter = new Parameter.Default(beta);
		cox = new CoxPartialLikelihood(parameter, data);
		System.err.println("C " + cox.getLogLikelihood());

		System.err.println();

	}

//	test <- read.table(header=T, sep = ",", text = "
//	start, length, event, x1, x2
//  0, 4,  1,0,0
//	0, 3,  1,2,0
//	0, 3,  0,0,1
//	0, 2,  1,0,1
//	0, 2,  1,1,1
//	0, 1,  0,1,0
//	0, 1,  1,1,0
//	")
//
//	gold <- coxph(Surv(length, event) ~ x1, test, ties = "breslow")
//	coef(gold)
//	logLik(gold)
//
//	data <- createCyclopsData(Surv(length, event) ~ x1, data = test,
//	modelType = "cox")
//
//	cyclops <- fitCyclopsModel(data)
//	coef(cyclops)
//	logLik(cyclops)
}
