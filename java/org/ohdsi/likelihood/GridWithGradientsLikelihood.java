/*******************************************************************************
 * Copyright 2024 Observational Health Data Sciences and Informatics
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

import org.ohdsi.data.GridWithGradientsData;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Likelihood;

/**
 * @author Martijn Schuemie
 */
public class GridWithGradientsLikelihood extends AbstractModelLikelihood {

	private static final long serialVersionUID = -934094651050122633L;
	private final Parameter beta;
	private final GridWithGradientsData data;
	private boolean likelihoodKnown;
	private double logLikelihood; 
	private boolean storedLikelihoodKnown;
	private double storedLogLikelihood;

	public GridWithGradientsLikelihood(Parameter beta, GridWithGradientsData data) {
		super("GridWithGradientsLikelihood");
		this.beta = beta;
		this.data = data;
		addVariable(beta);
		likelihoodKnown = false;
	}

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
		// Do nothing
	}

	@Override
	public Model getModel() {
		return this;
	}

	private double computeLikelihood() {
		int n = data.point.length;
		double x = beta.getParameterValue(0);

		// Extrapolate to the left (linear)
		if (x < data.point[0]) {
			return data.value[0] + (x - data.point[0]) * Math.max(data.derivative[0], 0);
		}

		// Extrapolate to the right (linear)
		if (x > data.point[n - 1]) {
			return data.value[n - 1] + (x - data.point[n - 1]) * Math.min(data.derivative[n - 1], 0);
		}

		// Find the interval containing x
		for (int i = 0; i < n - 1; i++) {
			if (x >= data.point[i] && x <= data.point[i + 1]) {
				double t = (x - data.point[i]) / (data.point[i + 1] - data.point[i]);

				double h00 = 2 * Math.pow(t, 3) - 3 * Math.pow(t, 2) + 1;
				double h10 = Math.pow(t, 3) - 2 * Math.pow(t, 2) + t;
				double h01 = -2 * Math.pow(t, 3) + 3 * Math.pow(t, 2);
				double h11 = Math.pow(t, 3) - Math.pow(t, 2);

				double dx = data.point[i + 1] - data.point[i];

				return h00 * data.value[i] +
						h10 * dx * data.derivative[i] +
						h01 * data.value[i + 1] +
						h11 * dx * data.derivative[i + 1];
			}
		}
		throw new IllegalArgumentException("x is out of range, and extrapolation failed.");
	}

	@Override
	public double getLogLikelihood() {
		if (!likelihoodKnown) {
			logLikelihood = computeLikelihood();
			likelihoodKnown = true;
		}
		return logLikelihood;
	}

	@Override
	public void makeDirty() {
		likelihoodKnown = false;
	}

	@Override
	protected void handleModelChangedEvent(Model arg0, Object arg1, int arg2) {
		// Do nothing
	}

	public static void main(String[] args) {
		double[] point = new double[] {1.1, 2.1};
		double[] value = new double[] {1, 1};
		double[] derivative = new double[] {0.1, -0.1};
		GridWithGradientsData data = new GridWithGradientsData(point, value, derivative);
		Parameter parameter = new Parameter.Default(0d);
		Likelihood likelihood = new GridWithGradientsLikelihood(parameter, data);


		double[] xs = new double[] {0, 1, 2, 3};
		for (double x: xs) {
			parameter.setParameterValue(0, x);
			likelihood.makeDirty();
			System.out.println(likelihood.getLogLikelihood());
		}

		// According to R: 0.890 0.990 1.009 0.910

	}
}
