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
 * @author Martijn J. Schuemie
 */
public class ConditionalPoissonLikelihood extends AbstractModelLikelihood {

	protected final Parameter beta;
	protected final Parameter nuisance;

	protected final SortedCoxData data;
	protected final int N;
	protected final int J;

//	private final double[] oddsLabeledExposedByStratum;

	protected final int[] countByStratum;

	private boolean likelihoodKnown;
	private boolean storedLikelihoodKnown;

	private double logLikelihood;
	private double storedLogLikelihood;

	public ConditionalPoissonLikelihood(Parameter beta, Parameter nuisance, SortedCoxData data) {
		this("conditionalPoissonLikelihood", beta, nuisance, data);
	}

	public ConditionalPoissonLikelihood(String name, Parameter beta, Parameter nuisance, SortedCoxData data) {

		super(name);

		this.beta = beta;
		this.nuisance = nuisance;

		this.data = data;
		this.N = data.strata[data.strata.length - 1];
		this.J = data.x.length / N;

		if (beta.getDimension() != J || data.x.length % N != 0) {
			throw new IllegalArgumentException("Illegal dimension for covariate matrix");
		}

		addVariable(beta);

		if (nuisance != null) {
			addVariable(nuisance);
		}

		checkStrata(data.strata, data.y.length);

//		oddsLabeledExposedByStratum = computeOddsLabeledExposed(data.x, data.strata);
		countByStratum = computeCount(data.y, data.strata);

		beta.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
				beta.getDimension()));
	}

	public Parameter getParameter() {
		return beta;
	}

	protected double calculateLogLikelihood() {

		final int[] strata = data.strata;
		final int[] y = data.y;
		final double[] x = data.x;
		final double[] weight = data.weight;

		final double oddsRatio = Math.exp(beta.getParameterValue(0));

		final double[] beta = this.beta.getParameterValues(); // TODO Faster as copy or direct access?

		final double sensitivity = nuisance.getParameterValue(0);
		final double reciprocalOddsSensitivity = (1.0 - sensitivity) / sensitivity;

		double logLikelihood = 0.0;

		int last = 0;
		for (int k = 0; k < strata.length; ++k) {

			double stratumDenominator = 0.0;
			for (int i = last; i < strata[k]; ++i) {

				double rowXBeta = 0.0;
				for (int j = 0; j < J; ++j) {
					rowXBeta += x[i * J + j] * beta[j]; // TODO Need to form xTranspose
				}

				logLikelihood += weight[i] * y[i] * rowXBeta;
				stratumDenominator += Math.exp(rowXBeta);
			}

			logLikelihood -= countByStratum[k] * Math.log(stratumDenominator);

			last = strata[k];
		}

		return logLikelihood;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object o, int i) { }

	@Override
	protected void handleVariableChangedEvent(Variable variable, int i, Variable.ChangeType changeType) {
		if (variable == beta) {
			likelihoodKnown = false;
		} else if (variable == nuisance) {
			likelihoodKnown = false;
		} else {
			throw new IllegalArgumentException("Unknown variable");
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

	private static int[] computeCount(int[] y, int[] strata) {
		int[] counts = new int[strata.length];

		int last = 0;
		for (int k = 0; k < strata.length; ++k) {
			int count = 0;
			for (int i = last; i < strata[k]; ++i) {
				count += y[i];
			}

			counts[k] = count;
			last = strata[k];
		}

		return counts;
	}

//	private static double[] computeOddsLabeledExposed(double[] x, int[] strata) {
//
//		double[] odds = new double[strata.length];
//
//		int last = 0;
//		for (int k = 0; k < strata.length; ++k) {
//
//			double sum = 0.0;
//			for(int i = last; i < strata[k]; ++i) {
//				sum += x[i];
//			}
//			double proportion = sum / (strata[k] - last);
//			odds[k] = proportion / (1.0 - proportion);
//
//			last = strata[k];
//		}
//
//		return odds;
//	}

	private static void checkStrata(int[] strata, int N) {

		boolean valid = true;

		for (int k = 1; k < strata.length; ++k) {
			if (strata[k] <= strata[k - 1]) {
				valid = false;
				break;
			}
		}

		if (strata[strata.length - 1] != N) {
			valid = false;
		}

		if (!valid) {
			throw new IllegalArgumentException("Invalid strata");
		}
	}

	public static SortedCoxData exampleData() {
		int[] y = { 0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,1,1,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,1,1,0,0,1,0,0,0,1,1,1,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,0,0,1,1,0,1,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,1,0,0,0,0,0,1,1,0,1,1,0,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,1,1,1,1,0,1,0,1,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,1,0,1,1,0,1,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,1,1,1,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,1,0,1,1,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,1,1,1,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,1,0,1,0,0,1,0,1,1,1,0,0,0,0,0,1,0,0,1,1,0,0,1,0,1,1,0,0,0,1,1,0,0,0,1,0,0,0,0,0,1,0,1,1,0,0,1,0,1,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1,0,0,1,0,0,0,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,0,1,1,0,1,0,1,1,0,1,1,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,1,0,1,0,0,0,0,0,1,0,1,1,0,1,0,0,1,0,0,0,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,0,1,1,0,0,0,0,1,0,0,1,1,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,1,0,1,1,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,1,0,1,1,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,1,1,0,1,1,0,0,0,1,1,1,0,1,1,0,0,1,1,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,1,1,0,0,1,1,1,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,1,0,1,1,0,0,0,0,0,1,0,0,1,0,1,1,0,1,0,0,1,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,0,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,1,0,1,1,0,0,0,0,0,1,1,0,1,0,1,0,1,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,0,1,0,1,1,1,1,0,1,1,1,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,1,1,0,0,0,0,1,0,0,0,1,0,1,0,1,1,0,1,0,1,1,1,1,1,0,1,0,1,0,0,1,0,1,0,1,0,0,1,0,0,0,0,1,1,1,0,1,0,0,1,0,1,1,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,0,1,0,0,0,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,1,1,1,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,0,0,1,0,1,1,0,0,0,0,0,0,1,1,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,1,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,1,0,0,0,0,0,1,0,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0,1,0,1,0,1,1 };
		double[] x = { 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,1,0,0,1,0,1,0,0,0,1,0,1,0,0,1,1,0,0,1,0,0,0,1,1,1,1,0,0,1,1,1,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,1,0,1,0,0,1,0,1,0,0,0,1,1,1,1,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,1,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,1,1,0,1,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,0,0,1,0,0,0,0,0,1,0,1,1,0,0,1,0,0,0,1,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,0,1,0,1,0,0,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,0,0,1,1,1,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,0,0,0,1,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,1,1,1,0,0,0,1,1,0,1,1,0,0,0,1,0,0,1,1,0,0,0,0,0,0,1,0,0,0,1,0,0,1,1,1,1,0,1,1,0,0,0,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,1,1,0,0,1,0,1,1,0,0,1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,0,1,1,0,1,1,0,0,0,1,1,1,1,0,1,1,1,1,0,0,1,1,0,1,0,0,1,0,0,1,1,1,0,1,1,0,1,0,1,0,1,1,0,1,1,1,0,1,0,0,0,1,0,1,1,0,1,0,0,0,0,0,1,1,0,1,1,1,1,0,1,1,0,0,1,1,1,0,0,1,0,0,1,1,1,0,0,1,1,1,0,1,0,0,0,0,0,0,1,0,1,1,0,0,1,0,0,1,1,1,1,1,0,1,0,1,1,0,0,1,1,1,0,1,0,1,1,1,1,0,0,0,1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,1,0,0,0,0,0,0,1,0,1,1,0,1,1,0,1,0,1,1,1,0,1,0,1,0,0,1,1,1,0,1,0,1,1,1,0,0,1,0,1,1,0,1,0,1,0,1,1,0,0,0,1,0,0,1,1,0,1,1,0,1,0,1,1,0,1,0,1,0,0,1,0,1,0,1,1,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,1,0,1,0,0,0,0,0,1,1,0,1,0,1,1,0,1,0,1,1,1,0,0,0,0,0,1,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,1,0,1,1,0,0,0,0,0,1,0,0,1,1,1,0,1,1,1,1,1,1,1,0,1,0,1,1,0,0,0,0,1,0,1,0,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,1,1,1,1,0,1,1,1,1,1,0,0,1,0,0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,1,0,0,1,1,1,0,1,0,0,1,0,1,0,0,1,1,0 };
		double[] weight = { 1.6,1,0.9,1.8,1.1,0.5,0.7,0.3,1.6,1.6,0.4,1.5,1.6,1,1.2,1,0.5,0.8,1.4,1,0.3,1.2,1.7,1.6,0.4,1,0.2,1.6,0.9,0.9,0.3,1.3,1.8,0.5,1.7,0.8,0.9,1.8,0.4,1.5,1.4,0.3,1.8,1,1,1.5,0.6,1,0.9,0.6,0.4,0.3,0.4,1.4,1.2,1.7,0.5,0.8,0.5,1.9,1.8,1,1.7,0.5,1.4,0.7,0.9,1.6,1.6,1.7,0.7,1.9,1.3,0.9,1.5,1.1,1,2,1,1.2,0.3,1.7,1.4,0.7,0.9,0.7,1.1,0.4,1.1,0.5,1.1,0.4,1.3,1.6,1.7,0.7,1.3,1.1,1.4,1.8,0.9,0.4,0.5,0.3,0.4,1.4,1.1,1.5,1.7,1.4,0.8,1.5,1.2,1.8,1.6,1,1.6,1,0.6,0.3,1.7,1,1,0.9,1.5,1.2,1.5,1.1,0.8,0.7,0.6,0.2,1.1,0.9,0.7,1.9,1.3,0.6,0.6,0.3,0.9,1,0.8,0.6,1.2,1.1,0.8,0.6,1.7,1.6,0.8,1.9,1.5,0.8,0.6,0.5,0.8,1.3,1.8,1.6,1.9,0.3,0.4,1.4,1.2,2,1.1,1.6,1.3,0.8,0.4,0.2,0.6,1.4,0.7,0.4,1,1.5,0.7,0.5,0.5,1.1,0.3,0.7,0.9,1.3,1.8,1.7,1.7,1.3,1.7,1.2,1.9,1.7,1.6,1.8,1.7,1,1.2,1.6,0.4,0.7,1.2,0.4,1.3,1.4,1.1,0.9,1,1.7,1.3,1.3,1.6,0.8,1.3,0.7,1.7,1.8,1.8,1.4,0.3,1.9,1.7,0.6,1.7,0.7,1.9,0.3,1.9,0.6,1.9,1.7,1.2,0.6,1.5,0.5,1.4,1.8,0.6,0.9,2,1.6,0.3,1.6,1.7,1.2,1.5,0.7,1.9,0.8,1.6,1.2,0.6,0.9,1,1.2,1.9,1.6,1.8,1.3,1.8,0.9,1.1,1.5,0.7,1.8,1.1,2,0.6,1,1.8,1.7,0.3,0.8,1.6,1.2,1.6,0.7,1.3,0.4,0.9,1.7,0.5,1.9,0.2,0.2,0.4,0.4,0.9,1.3,0.9,0.3,1.2,0.7,0.9,1.1,1.8,1.3,1.1,1.1,2,0.9,1.4,0.9,1.4,2,0.8,0.5,0.9,0.6,1.1,1.3,1.8,1.4,1.9,1.6,1.4,0.8,1,0.5,1.1,0.3,1.2,1.1,1.3,0.5,0.4,0.4,0.6,1.5,0.4,1.7,1.7,0.7,1.1,1.9,1,1.7,1.2,1.1,1,0.6,1.7,1,1.3,0.7,1.7,0.9,0.3,1.2,1.8,0.6,1.1,1.4,1.1,0.3,0.3,0.3,1.2,1.7,1.3,0.6,1,0.2,1.2,0.7,1.2,0.5,0.9,0.2,1.2,0.4,1.6,1.2,0.3,0.2,0.5,1.3,0.4,1.3,0.8,1.4,1.8,0.7,1,1.4,1.5,1.6,0.5,1.4,0.9,1.1,0.3,1.6,1.7,1.7,0.8,1.6,0.6,0.4,1.6,1,0.3,1.8,0.6,0.3,1.6,1.6,1.5,2,1.1,1.3,1.5,1.4,1.8,0.5,1.1,1.3,1.7,1.5,1.8,0.7,1.9,0.9,0.5,0.6,0.5,1.5,1,0.2,1.2,0.9,0.3,0.9,0.6,0.9,0.3,0.4,1.8,1.7,0.3,1.9,1.9,0.9,0.5,0.5,1.7,0.3,1.1,1.3,1.3,1.1,1.4,1.5,1.1,0.8,0.8,1,1.6,0.9,1.7,1.2,0.2,1.7,0.5,0.6,0.8,1.7,1.1,1.5,0.9,0.8,1.6,0.6,0.6,1.7,0.8,1.5,0.3,0.5,1.7,0.8,1.4,1.7,1.8,1,1.4,1.2,0.9,0.4,1.2,1.9,0.5,1.3,1.1,1.4,1.8,1.7,2,1.8,1,0.4,0.9,1.7,1.7,0.8,1.3,1.6,1.6,0.9,0.4,0.5,1.9,0.2,1.8,1.9,1.3,1.1,0.7,1,0.7,1.6,1.3,1.5,1.2,0.4,0.7,1.6,1,0.9,1.6,0.9,0.3,1.9,0.6,0.6,1.1,0.4,1.2,1.2,1.5,1.6,1.6,0.5,1,0.7,1.4,1.5,0.3,0.7,0.3,1.2,0.4,1.9,1.8,1.7,1.7,1,0.9,1.8,0.9,1.5,1.2,1.2,0.7,0.5,1.6,1.1,1.3,1.1,0.5,0.9,1.9,1.2,0.6,0.4,1,1.9,1.9,1.4,1.6,1.8,1.8,0.8,0.6,1.1,0.7,1.6,0.4,1.1,0.7,2,1.1,0.5,1.1,1.2,0.7,0.7,0.4,1.3,0.7,1.3,0.4,0.7,1.7,0.9,1,0.3,1.1,1.1,1.4,1.1,1.8,1.7,1.3,0.5,1.3,1.1,0.8,1.8,1.4,1.1,0.5,1.8,1,1.6,0.9,1.9,1.5,0.9,1.8,1.3,1.1,1.3,0.4,0.3,0.8,0.9,1.7,0.5,0.3,1.1,1.5,1.4,1.7,0.4,0.3,1.5,0.5,0.7,1.2,2,0.8,1.6,1.3,1.2,0.7,0.3,1.1,0.7,1.5,1.3,1.2,1.2,1.7,1,0.2,1.5,0.9,1.3,1.9,0.9,1,1.7,1.5,0.4,0.3,1.3,1.4,1.8,2,1.4,1.2,1.8,1.6,1,1.2,1.5,0.4,0.4,1,1.2,1.5,1,1.4,0.9,1,0.4,0.3,0.2,1.4,1.4,1.5,1.8,1.6,0.9,1.9,1.6,1.1,1.4,0.5,1.4,0.9,0.8,1.8,0.4,0.3,0.7,1.2,1.3,0.9,2,1.7,1.3,0.9,0.8,0.6,0.4,0.9,1.2,1.3,0.5,0.6,0.4,1.1,1.1,1.7,0.6,1.4,2,0.4,0.2,1.7,1.2,2,1.1,1.2,1,0.5,1.9,0.5,0.3,1.2,1.2,2,0.3,1.2,1.9,0.5,0.8,1.3,0.9,1.2,1,0.9,1.8,0.2,0.2,1.1,0.5,1,1.8,1.4,0.5,0.3,0.6,1.8,1.6,1.7,0.4,0.3,0.5,1.5,0.7,1.5,0.6,1.6,0.5,1.3,1.3,1.9,1.5,1.4,1.1,1.2,2,1.4,2,0.9,1.9,1.1,0.9,2,1.6,1.4,0.5,1.3,1.6,0.5,1.9,2,1.2,1.1,1.6,0.5,1.8,2,0.7,0.4,1.7,0.9,1.8,0.4,1.1,1.7,0.5,1.5,2,1,1.7,0.7,1.2,1.6,0.4,0.6,0.4,1.3,1.9,1.4,1,1.1,1,1.9,1.3,0.7,0.6,1,1,1,1.5,1.6,0.9,1.9,1.1,0.4,0.4,1.1,1.9,1,0.2,1.2,0.5,1,1.8,1.4,1,1.8,0.5,1.5,1.9,0.6,0.8,0.2,0.9,0.8,0.3,1,0.3,0.4,1.4,1.4,1,1.1,0.9,1,0.7,0.6,1.6,1.5,1,0.4,0.8,0.8,1,1.2,1.1,0.5,1.6,1.1,1.8,2,1.5,0.3,1,1.4,1.2,1.4,0.9,1,1.6,1,1.5,0.5,1.1,1.2,1.4,1,1.8,0.7,1.9,1.8,0.4,0.3,0.7,0.7,0.4,0.7,1.4,0.2,1.5,0.5,1.4,1.7,1.7,0.8,0.3,0.6,0.8,0.6,0.5,0.8,1.1,0.7,1.9,0.8,0.4,0.4,0.7,0.4,0.7,0.5,0.6,1.2,1.1,0.7,0.3,0.6,1,0.6,1.9,0.5,0.5,0.8,0.7,1.9,0.7,1.2,0.3,1.4,0.2,1,0.9,1.5,1.7,1.8,0.7,1.7,1.6,2,1.1,0.5,0.5,1.2,0.4,1.8,0.3,0.5,0.2,1.2,0.6,1.9,0.6,0.6,1.8,0.5,0.3,0.3,1,0.8,1.1,1.4,1.8,0.2,1.7,1.9,0.6,0.3,0.4,1.6,1.6,1.6,0.9,0.3,1.7,0.2,0.5,0.3,1.7,1.8,0.6,1.6,1.1,1.1,0.9,1.3,0.5,1.9,0.2,1.5,1.2,1.5,1.7,1.5,1.6,2,1,1.2,0.8,0.8,1,1.1,0.4,0.8,1.9,1.6,1.7,1.4,0.5,0.8,1.4,1.3,1.4,0.2,0.4,0.8,0.9,1.5,1,1.2,0.3,1.1,0.6,0.4,0.9,1.6,0.9,0.4,1,1.6,1.7,1.9,1.2,0.3,1.8,0.8,1.5,1,1.6,0.5,1.5,1.6,1.6,1.1,1,1.1,0.2,0.4,1.3,0.4,1.8,0.4,1.2,1.7,1.7,0.7,1.3,1.9,2,0.4,0.4,0.3,1.4,1,1.9,0.5,0.5,0.7,0.3,1.1,0.4,0.8,1.4,0.9,1.5,1.6,1.6,1.2,1.6,1.4,1.4,0.8,1.9,1.7,0.8,0.8,1.6,1,1.1,1.4,1.1,0.3,0.7,1.8,2,1.1,1.8,1.2,0.6,1.8,1.9,0.9,0.2,1.3,0.9,1.1,0.9,0.4,1.9,0.6,0.8,1.1,0.4,0.8,0.9,0.9,1.8,1.6,1,0.9,0.9,1.8,1.7,0.9,1.5,1.2,1.3,0.3,1.5,1.4,1.7,0.4,0.7,0.3,1.2,1.8,1.1,0.3,0.4,0.5,0.4,1.9,0.5,0.5,0.4,1.8,1.2,0.5,0.9,1.9,0.9,1.7,0.5,1.4,1.1,1.2,0.9,1.7,0.3,1.6,1.3,1.4,0.5,0.3,1.5,1.6,1.3,1.5,1,1.7,1.7,0.4,0.4,1.7,1.8,0.9,1.4,1.9,1.8,1,1.1,0.8,2,1.8,1,1.9,0.8,0.2,1.1,0.8,1.6,1.5,0.6,1.5,1,0.9,1.3,0.8,0.3,1.2,1.1,0.7,1,1,1.7,0.5,2,1.9,0.5,0.5,0.9,0.6,1.9,0.3,0.9,0.5,1.4,1.3,0.3,0.3,1.6,0.6,1.7,1.1,0.8,0.7,1.9,0.7,1.5,1.4,1.9,1.6,1.3,1.5,0.8,1.8,0.6,1.7,0.6,1.4,1.7,1.8,0.2,0.6,1.5,1.6,0.9,1.7,0.9,1.7,1.9,0.8,1.1,0.4,1.9,1.5,1.7,1.4,0.6,1.6,1.6,1.7,1.6,0.4,0.5,1.8,1.9,0.3,0.2,0.6,0.8,1.2,0.6,0.2,1,1.3,1.8,1.6,1.2,1.1,0.5,1.9,0.5,1.6,0.9,0.5,0.4,1,0.5,0.3,1.7,0.4,0.4,1.7,1.2,2,1.9,1.8,1.4,2,1.1,2,0.9,0.6,0.7,0.6,1.2,0.7,1.8,1.9,0.7,0.6,1.2,1.4,1.6,1.6,2,0.3,1.4,1.1,1.7,1.9,0.6,1.5,2,0.7,1.2,1.5,1.2,0.5,0.6,1.8,1,0.4,1,1.3,0.8,0.4,1,1.1,1.4,1.1,0.7,1.5,1.2,0.4,1.2,1,1.9,1.1,1.1,2,1.4,1.7,1.7,1.2,0.2,0.9,1.9,1.7,1.8,0.6,1.4,1,0.9,1.4,1.7,0.4,1.8,0.8,0.9,0.4,0.9,1.7,1.5,0.8,0.7,0.4,1.9,0.4,1.6,0.7,0.5,1.2,1.3,0.7,1,1.3,0.4,1.6,0.5,0.9,0.4,1.9,0.6,0.8,0.8,0.7,1.6,1.1,1.2,0.7,0.9,1.9,0.9,1.9,0.8,1.5,0.4,0.2,1.6,0.6,0.6,1,1.7,1.3,0.9,1,1.8,1.6,0.2,1,1.7,0.9,0.6,0.6,0.5,0.4,0.9,0.9,0.5,0.9,0.4,0.9,1,1.1,0.8,0.8,1.2,0.5,1.6,1.7,1.6,0.6,0.6,0.6,0.6,0.5,1.5,0.9,0.4,1.1,1.3,0.6,1,1.6,0.8,1.2,1.1,0.7,0.4,0.7,0.8,1.9,0.2,1.6,0.8,2,0.8,1.8,0.9,1.6,0.9,0.4,1.4,0.6,1,0.7,1.5,1.1,0.6,1.8,1.4,1,1.3,1.2,0.2,1.2,0.6,1.6,1.7,1.6,1.3,0.7,1.9,0.5,1,0.3,2,2,0.7,0.3,0.4,1.5,0.7,0.2,0.9,0.3,0.6,1,2,0.4,0.3,1.2,2,2,1.7,1,0.3,1.4,0.3,0.8,0.5,1.5,1.8,0.6,1.2,0.5,1.2,1,1.7,1.6,1.1,1.4,1.5,0.7,0.6,1,0.9,1.8,1.9,0.4,1.6,1.2,1,1.4,1.7,1.5,1.3,1,1.4,1.5,0.5,1.2,1.7,1.2,0.4,0.4,0.3,1.7,1.6,0.3,1.5,1.9,1.8,0.2,1.3,1.1,0.8,1.8,1.3,0.7,1,1.5,0.9,1.9,1,2,1.2,1.6,1.4,1.7,0.5,0.8,0.7,1.3,1.4,2,1.7,1.3,1.1,0.3,0.8,0.3,0.4,1.1,1.4,0.4,1.3,0.5,0.7,0.4,2,1,1.8,1.5,0.6,1.1,0.5,0.6,0.2,0.9,1.4,0.6,0.2,0.9,1.5,0.4,2,0.3,0.7,0.8,1.8,1,0.8,1.6,0.8,0.8,1.2,1,1.1,1.2,0.4,0.4,0.6,0.7,1.3,1.6,1.1,0.4,1.4,1.2,1.4,1.9,0.4,1.6,0.8,1.6,1.5,1.6,2,1.8,1.8,1.5,1.6,1.5,0.5,0.9,1.4,0.4,1.9,1.2,1.5,1.9,0.6,0.8,1,1.3,0.6,1.9,0.4,0.5,1.4,1.3,0.6,1.7,1.6,0.3,1.9,0.6,1.3,0.6,1.4,0.3,1.9,0.6,0.3,1.4,0.8,1.9,1.6,1.2,1.1,0.3,1.6,1.8,0.3,1.7,1.6,1.6,0.5,0.5,1.1,0.3,0.7,1.5,1.8,0.4,0.6,0.8,0.4,2,1.2,1,1.4,1,1.2,1,1.6,0.7,0.7,1.5,0.6,1.8,1.9,0.4,1.6,1.1,1.2,1.1,1,0.6,0.6,0.3,1.5,0.3,0.7,0.5,0.3,1.8,1.1,1.8,1.8,1.8,0.5,1.6,1.9,0.5,1.1,1.6,0.5,0.2,1,0.8,0.3,0.6,0.3,0.4,0.2,0.6,1,1.2,0.8,0.9,1,0.7,1.7,1,0.7,0.6,0.2,1.3,0.6,0.7,0.7,1.3,0.3,0.3,1.9,0.4,1.1,0.5,0.4,0.6,1.4,1.1,1.4,2,0.4,0.4,1.9,0.5,1.2,0.3,1.6,0.4,0.3,1.1,0.8,0.9,1.5,0.6,0.8,0.9,1.4,1.6,0.4,0.5,1.3,1.9,0.7,1.7,1.4,1.5,1,1.2,1.3,1.7,0.3,0.3,0.4,1.8,2,2,0.3,0.6,0.2,0.4,0.7,1,1.2,0.7,1.7,0.9,1.8,1.2,0.3,1.1,1.6,1.4,0.7,1.1,0.3,0.4,1.9,1,1,1.7,1.8,1.7,1.1,1.6,1.9,0.8,1.1,1.8,0.9,0.3,1.4,0.9,0.4,1.4,1.3,1.3,0.5,1.2,1.1,1.6,1.4,1,1.7,1.6,1.1,0.5,1.6,0.9,0.7,1.5,0.8,1.9,0.4,0.7,1.4,1.8,1.1,2,0.4,0.6,2,1.9,1.2,0.6,0.7,1,1.4,1.2,2,0.4,1.8,1,1,0.7,0.5,0.8,0.4,1.6,0.6,1.5,1.4,0.5,0.8,0.4,1.9,1.9,0.6,1.6,0.2,1.8,1.8,1.1,0.5,1.1,1.8,1.9,2,0.9,0.9,1.2,1.6,1.4,1.4,1.6,0.3,1.9,0.6,0.4,0.8,0.3,1.6,1.4,0.9,1.4,1.1,0.5,0.4,1.8,0.6,1.8,1.8,1.2,1.5,0.8,0.3,1.8,1,1.5,1.4,1.3,1.1,0.7,0.4 };
		int[] strata = { 400,801,1207,1601,2000 };

		return new SortedCoxData(y, x, strata, null, weight);
	}

	public static void main(String[] args) {

		double beta = 0.1;
		double sensitivity = 0.9;

		SortedCoxData data = exampleData();

		Parameter parameter = new Parameter.Default(new double[] { beta });

		Parameter nuisance = new Parameter.Default(new double[] { sensitivity });

		Likelihood cox = new ConditionalPoissonLikelihood(parameter, nuisance, data);

		System.err.println(cox.getLogLikelihood());
	}
}
