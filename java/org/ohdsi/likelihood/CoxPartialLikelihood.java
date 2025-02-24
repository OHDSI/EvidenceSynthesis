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

import org.apache.commons.math.util.FastMath;
import org.ohdsi.data.SortedCoxData;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc A. Suchard
 */
public class CoxPartialLikelihood extends AbstractModelLikelihood {

	private static final long serialVersionUID = 5671381056989317642L;
	private final Parameter beta;
	private final SortedCoxData data;
	private final int N;
	private final double numeratorConstant;

	public CoxPartialLikelihood(Parameter beta, SortedCoxData data) {

		super("coxPartialLikelihood");
		this.beta = beta;
		this.data = data;
		this.N = data.y.length;

		addVariable(beta);

		final int[] y = data.y;
		final double[] x = data.x;

		double sum = 0.0;
		for (int i = 0; i < N; ++i) {
			sum += (y[i] == 1) ? x[i] : 0.0;
		}
		numeratorConstant = sum;
	}

	private double numeratorContribution() {
		return numeratorConstant * beta.getParameterValue(0);
	}

	private double denominatorContribution() {

		final double[] weight = data.weight;
		final double[] x = data.x;
		final int[] strata = data.strata;
		final double b = beta.getParameterValue(0);

		double sum = 0.0;

		int resetIndex = 0;
		double denominatorTotal = 0.0;
		for (int i = 0; i < N; ++i) {
			if (i == strata[resetIndex]) {
				denominatorTotal = 0.0;
				++resetIndex;
			}
			denominatorTotal += FastMath.exp(x[i] * b);
			double increment =  weight[i] * FastMath.log(denominatorTotal);
			sum += increment;

			//            System.err.println("\t\t" + i + " " + denominatorTotal + ":" + increment);
		}

		return sum;
	}

	private double calculateLogLikelihood() {
		return numeratorContribution() - denominatorContribution();
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object o, int i) { }

	@SuppressWarnings("rawtypes")
	@Override
	protected void handleVariableChangedEvent(Variable variable, int i, Variable.ChangeType changeType) {
		if (variable == beta) {
			likelihoodKnown = false;
		}
	}

	@Override
	protected void storeState() {
		storedLikelihoodKnown = likelihoodKnown;
		storedLogLikelihood = logLikelihood;
	}

	@Override
	protected void restoreState() {
		likelihoodKnown = storedLikelihoodKnown;
		logLikelihood = storedLogLikelihood;
	}

	@Override
	protected void acceptState() {
	}

	@Override
	public Model getModel() {
		return this;
	}

	@Override
	public double getLogLikelihood() {
		if (!likelihoodKnown) {
			logLikelihood = calculateLogLikelihood();
			likelihoodKnown = true;
		}

		return logLikelihood;
	}

	@Override
	public void makeDirty() {
		likelihoodKnown = false;
	}

	private boolean likelihoodKnown;
	private boolean storedLikelihoodKnown;

	private double logLikelihood;
	private double storedLogLikelihood;

	public static void main(String[] args) {

		// No strata
		//   		int[] y = new int[] { 1, 1, 0, 1, 1, 0, 1 };
		//   		double[] x = new double[] { 0, 2, 0, 0, 1, 1, 1 };
		//   		int[] strata = new int[] { 7 };
		//   		double beta = 0.3883064;

		// With strata
		int[] y = new int[] { 1, 1, 0, 1, 0, 1, 1 };
		double[] weight = new double[] { 1, 1, 0, 1, 0, 1, 1 };
		double[] x = new double[] { 0, 2, 1, 1, 0, 0, 1 };
		int[] strata = new int[] { 4, 7 };
		double beta = 1.205852;
		// logLike = -2.978028

		SortedCoxData data = new SortedCoxData(y, x, strata, weight);

		Parameter parameter = new Parameter.Default(beta);

		Likelihood cox = new CoxPartialLikelihood(parameter, data);

		System.err.println(cox.getLogLikelihood());
	}
}
