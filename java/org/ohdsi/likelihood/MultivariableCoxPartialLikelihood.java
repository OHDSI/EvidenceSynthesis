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
import org.apache.commons.math.util.FastMath;
import org.ohdsi.data.ColumnMajorSortedCoxData;
import org.ohdsi.data.CoxData;
import org.ohdsi.data.SortedCoxData;

/**
 * @author Marc A. Suchard
 */
public class MultivariableCoxPartialLikelihood extends ConditionalPoissonLikelihood {

	// Currently assumes that data.x is row-major

	public MultivariableCoxPartialLikelihood(Parameter beta, SortedCoxData data) {
		super("multivariableCoxPartialLikelihood", beta, null, data);
	}

	@Override
	protected double calculateLogLikelihood() {

		final int[] y = data.y;
		final double[] x = data.x;
		final int[] n = data.n;
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
				rowXBeta += x[i * J + j] * beta[j];
			}

			denominator += FastMath.exp(rowXBeta); // TODO Implement weights
			logLikelihood +=  y[i] * rowXBeta - n[i] * FastMath.log(denominator); // TODO Could pre-compute total numerator as function of beta
		}

		return logLikelihood;
	}

	public static SortedCoxData exampleMultivariableData() {
		int[] y = new int[] { 1, 1, 0, 1, 1, 0, 1 };
		double[] x = new double[] { 0, 2, 0, 0, 1, 1, 1,
									0, 0, 1, 1, 1, 0, 0 };
		int[] yTimesTieCount = new int[] { 1, 1, 0, 1, 1, 0, 1 };
		int[] strata = new int[] { 7 };
		return new ColumnMajorSortedCoxData(y, x, strata, yTimesTieCount, null);
	}

	public static SortedCoxData exampleBladder() {
		int[] y = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0 };
		double[] time = { 1,1,1,1,4,4,4,4,7,7,7,7,10,10,10,10,6,10,10,10,14,14,14,14,18,18,18,18,5,18,18,18,12,16,18,18,23,23,23,23,10,15,23,23,3,16,23,23,3,9,21,23,7,10,16,24,3,15,25,25,26,26,26,26,1,26,26,26,2,26,26,26,25,28,28,28,29,29,29,29,29,29,29,29,29,29,29,29,28,30,30,30,2,17,22,30,3,6,8,12,12,15,24,31,32,32,32,32,34,34,34,34,36,36,36,36,29,36,36,36,37,37,37,37,9,17,22,24,16,19,23,29,41,41,41,41,3,43,43,43,6,43,43,43,3,6,9,44,9,11,20,26,18,48,48,48,49,49,49,49,35,51,51,51,17,53,53,53,3,15,46,51,59,59,59,59,2,15,24,30,5,14,19,27,2,8,12,13,1,1,1,1,1,1,1,1,5,5,5,5,9,9,9,9,10,10,10,10,13,13,13,13,3,14,14,14,1,3,5,7,18,18,18,18,17,18,18,18,2,19,19,19,17,19,21,21,22,22,22,22,25,25,25,25,25,25,25,25,25,25,25,25,6,12,13,26,6,27,27,27,2,29,29,29,26,35,36,36,38,38,38,38,22,23,27,32,4,16,23,27,24,26,29,40,41,41,41,41,41,41,41,41,1,27,43,43,44,44,44,44,2,20,23,27,45,45,45,45,2,46,46,46,46,46,46,46,49,49,49,49,50,50,50,50,4,24,47,50,54,54,54,54,38,54,54,54,59,59,59,59 };
		double[] x = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,1,1,1,1,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1,3,3,3,3,1,1,1,1,2,2,2,2,1,1,1,1,4,4,4,4,2,2,2,2,4,4,4,4,2,2,2,2,1,1,1,1,6,6,6,6,5,5,5,5,1,1,1,1,3,3,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,1,1,1,1,1,1,1,1,2,2,2,2,1,1,1,1,6,6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,1,1,1,1,7,7,7,7,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1,2,2,2,2,1,1,1,1,1,1,1,1,6,6,6,6,3,3,3,3,1,1,1,1,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,5,5,5,5,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,4,4,4,4,4,4,4,4,3,3,3,3,1,1,1,1,1,1,1,1,4,4,4,4,1,1,1,1,3,3,3,3 };

		return new CoxData(y, time, x).getSortedData();
	}

	public static void main(String[] args) {

		int[] y;
		double[] x;
		int[] yTimesTieCount;
		int[] strata;
		double beta;
		double[] betas;

		SortedCoxData data;
		Parameter parameter;
		Likelihood cox;

		// No ties, no strata, no weights
		y = new int[] { 1, 1, 0, 1, 1, 0, 1 };
		x = new double[] { 0, 2, 0, 0, 1, 1, 1 };
		yTimesTieCount = new int[] { 1, 1, 0, 1, 1, 0, 1 };
		strata = new int[] { 7 };
		beta = 0.3883064;
		// logLik = -5.401371

		data = new SortedCoxData(y, x, strata, yTimesTieCount, null);
		parameter = new Parameter.Default(beta);
		cox = new MultivariableCoxPartialLikelihood(parameter, data);
		System.err.println("M " + cox.getLogLikelihood());
		data = new SortedCoxData(y, x, strata, yTimesTieCount, null);
		parameter = new Parameter.Default(beta);
		cox = new CoxPartialLikelihood(parameter, data);
		System.err.println("C " + cox.getLogLikelihood());

		System.err.println();

		// No ties, with strata, no weights
		y = new int[] { 1, 1, 0, 1, 0, 1, 1 };
		yTimesTieCount = new int[] { 1, 1, 0, 1, 0, 1, 1 };
		x = new double[] { 0, 2, 1, 1, 0, 0, 1 };
		strata = new int[] { 4, 7 };
		beta = 1.205852;
		// logLike = -2.978028

		data = new SortedCoxData(y, x, strata, yTimesTieCount, null);
		parameter = new Parameter.Default(beta);
		cox = new MultivariableCoxPartialLikelihood(parameter, data);
		System.err.println("M " + cox.getLogLikelihood());
		data = new SortedCoxData(y, x, strata, yTimesTieCount, null);
		parameter = new Parameter.Default(beta);
		cox = new CoxPartialLikelihood(parameter, data);
		System.err.println("C " + cox.getLogLikelihood());

		System.err.println();

		// No ties, no strata, no weights, 2 covariates
		data = exampleMultivariableData();
		betas = new double[] { 0.823619, 1.518213 };
		// logLike = -4.879409

		parameter = new Parameter.Default(betas);
		cox = new MultivariableCoxPartialLikelihood(parameter, data);
		System.err.println("M " + cox.getLogLikelihood());

		// Larger `survival::bladder` example
		cox = new MultivariableCoxPartialLikelihood(
				new Parameter.Default(new double[] { -0.4608773, -0.1012988 }),
				exampleBladder());
		System.err.println("M " + cox.getLogLikelihood());
		// logLike = -596.3328
	}

//	test <- read.table(header=T, sep = ",", text = "
//	start, length, event, x1, x2
//	0, 4,  1,0,0
//	0, 3.5,1,2,0
//	0, 3,  0,0,1
//	0, 2.5,1,0,1
//	0, 2,  1,1,1
//	0, 1.5,0,1,0
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
