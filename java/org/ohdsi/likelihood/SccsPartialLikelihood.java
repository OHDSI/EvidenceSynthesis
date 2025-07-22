/*******************************************************************************
 * Copyright 2025 Observational Health Data Sciences and Informatics
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

import java.util.Arrays;

import org.apache.commons.math.util.FastMath;
import org.ohdsi.data.SccsData;

import com.github.lbfgs4j.LbfgsMinimizer;
import com.github.lbfgs4j.liblbfgs.Function;
import com.github.lbfgs4j.liblbfgs.Lbfgs;
import com.github.lbfgs4j.liblbfgs.LbfgsConstant.LBFGS_Param;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Martijn Schuemie
 */
public class SccsPartialLikelihood extends AbstractModelLikelihood {

	private static final long serialVersionUID = 5911070778889767445L;
	private final Parameter beta;
	private final SccsData data;
	private final double n;
	private final double[] exps;
	private final double[] xps;
	private final int[] idxs;
	private boolean likelihoodKnown;
	private double minLogLikelihood; // Internally using min log likelihood because lbfgs4j minimizes
	private boolean storedLikelihoodKnown;
	private double storedMinLogLikelihood;
	private double[] minGradient; // Internally using min gradient because lbfgs4j minimizes

	public SccsPartialLikelihood(Parameter beta, SccsData data) {

		super("SccsPartialLikelihood");
		this.beta = beta;
		this.data = data;
		this.n = data.y.length;

		addVariable(beta);

		int stratumSize = 1;
		int maxStratumSize = 0;
		for (int i = 1; i < n; i++) {
			if (data.stratumId[i-1] != data.stratumId[i]) {
				if (stratumSize > maxStratumSize)
					maxStratumSize = stratumSize;
				stratumSize = 0;
			}
			if (data.y[i] != 0)
				stratumSize++;

		}
		if (stratumSize > maxStratumSize)
			maxStratumSize = stratumSize;

		exps = new double[maxStratumSize];
		xps = new double[maxStratumSize];
		idxs = new int[maxStratumSize];
		minGradient = new double[data.x[0].length];
		likelihoodKnown = false;
	}

	private void computeLogLikelihoodAndGradient(double pA, double[] pX) {
		// Martijn sucks at math, so used this Wolfram Alpha query to get formula for gradient:
		// differentiate y*log(Subscript[v,1]*e^Subscript[x,11]*a + Subscript[x,12]*b/Subscript[v,1]*e^Subscript[x,11]*a + Subscript[x,12]*b + Subscript[v,2]*e^Subscript[x,21]*a + Subscript[x,22]*b) wrt b

		int cursor = 0; 
		double sumExps = 0;
		double[] sumExpXs = new double[pX.length];
		Arrays.fill(sumExpXs, 0);
		minLogLikelihood = 0;
		Arrays.fill(minGradient, 0);
		for (int i = 0; i < n; i++) {
			double xp = data.a[i] * pA;
			for (int j = 0; j < pX.length; j++)
				xp += data.x[i][j] * pX[j];
			double exp = data.time[i] * FastMath.exp(xp);
			sumExps += exp;
			for (int j = 0; j < pX.length; j++)
				sumExpXs[j] += data.x[i][j] * exp;
			if (data.y[i] != 0) {
				xps[cursor] = xp;
				exps[cursor] = exp;
				idxs[cursor] = i;
				cursor++;
			}			
			if (i == n-1 || data.stratumId[i] != data.stratumId[i+1]) {
				// End of stratum
				double logDenominator = FastMath.log(sumExps);
				for (int j = 0; j < cursor; j++) {
					int idx = idxs[j];
					minLogLikelihood -= data.y[idx]*(FastMath.log(exps[j])-logDenominator);
					for (int k = 0; k < pX.length; k++) {
						double part1 = data.y[idx] * FastMath.exp(-xps[j]) * sumExps / data.time[idx];
						double part2 = ((data.x[idx][k] * exps[j]) / sumExps); 
						double part3 = ((exps[j] * sumExpXs[k])/FastMath.pow(sumExps, 2));
						minGradient[k] -= part1 * (part2 - part3);
					}
				}
				cursor = 0;
				sumExps = 0;
				Arrays.fill(sumExpXs, 0);
			}
		}
	}

	private void refitModelAtNewBeta() {
		OptimizableFuntion function = new OptimizableFuntion();
		LBFGS_Param params = Lbfgs.defaultParams(); 
		params.epsilon = 1e-06;
		params.delta = 1e-06;
		LbfgsMinimizer minimizer = new LbfgsMinimizer(params, false);
		minimizer.minimize(function); 
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
		storedMinLogLikelihood = minLogLikelihood;
	}

	@Override
	protected void restoreState() {
		likelihoodKnown = storedLikelihoodKnown;
		minLogLikelihood = storedMinLogLikelihood;
	}

	@Override
	protected void acceptState() {
		// Do nothing
	}

	@Override
	public Model getModel() {
		return this;
	}

	@Override
	public double getLogLikelihood() {
		if (!likelihoodKnown) {
			refitModelAtNewBeta();
			likelihoodKnown = true;
		}
		return -minLogLikelihood;
	}

	@Override
	public void makeDirty() {
		likelihoodKnown = false;
	}

	@Override
	protected void handleModelChangedEvent(Model arg0, Object arg1, int arg2) {
		// Do nothing
	}

	class OptimizableFuntion implements Function {

		double[] lastPoint;

		public OptimizableFuntion() {
			lastPoint = new double[minGradient.length];
			Arrays.fill(lastPoint, 9999);
		}

		private void updateIfNeeded(double[] point) {
			boolean needUpdate = false;
			for (int i = 0; i < point.length; i++)
				if (lastPoint[i] != point[i]) {
					needUpdate = true;
					break;
				}
			if (needUpdate) {
				computeLogLikelihoodAndGradient(beta.getParameterValue(0), point);
			}				
		}

		@Override
		public int getDimension() {
			return minGradient.length;
		}

		@Override
		public double[] gradientAt(double[] arg0) {
			updateIfNeeded(arg0);
			return minGradient;
		}

		@Override
		public double valueAt(double[] arg0) {
			updateIfNeeded(arg0);
			return minLogLikelihood;
		}
	}
	
	public static void main(String[] args) {
		// Check gradient:
		//		double EPS = 0.00001;
		//		int[] y = {1, 0, 0, 1};
		//		double[] a = {1, 0, 1, 0};
		//		double[][] x = {{0.1, 0.2}, {0.2, 0.1}, {0.3, 0.1}, {0.5, 0.3}};
		//		double[] time = {3, 5, 3, 1};
		//		
		//		SccsData sccsData = new SccsData(y, a, x, stratumId, time);
		//		
		//		SccsPartialLikelihood spl = new SccsPartialLikelihood(null, sccsData);
		//		spl.computeLogLikelihoodAndGradient(0.1, new double[]{0.2, 0.3});
		//		double grad0 = spl.gradient[0];
		//		double grad1 = spl.gradient[1];
		//		double ll = spl.logLikelihood;
		//		System.out.println(ll);
		//		spl.computeLogLikelihoodAndGradient(0.1, new double[]{0.2+EPS, 0.3});
		//		double grad0Dm = (spl.logLikelihood - ll) / EPS;				
		//		System.out.println(grad0);
		//		System.out.println(grad0Dm);
		//		
		//		
		//		spl.computeLogLikelihoodAndGradient(0.1, new double[]{0.2, 0.3+EPS});
		//		double grad1Dm = (spl.logLikelihood - ll) / EPS;				
		//		System.out.println(grad1);
		//		System.out.println(grad1Dm);

		// Oxford example from Farrington and Whitaker (taken from Cyclops' unit tests)
		int[] y = new int[] { 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0 };
		double[] time = new double[] { 107, 21, 54, 183, 41, 21, 120, 183, 78, 21, 83, 183, 82, 21, 79, 183, 81, 21, 80, 183, 44, 21, 117, 183, 119, 21, 42, 183, 145, 21, 16, 183, 77, 21, 84, 183, 182, 183 };
		double[] a = new double[] { 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0 };
		double[][] x = new double[][] { {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {1}};
		int[] stratumId = new int[] { 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10 };
		SccsData data = new SccsData(y, a, x, stratumId, time);
		
		Parameter parameter = new Parameter.Default(-9999.9);// Some extreme value
		Likelihood sccs = new SccsPartialLikelihood(parameter, data);
		System.err.println(sccs.getLogLikelihood());
		// NaN
		
		parameter.setParameterValue(0, 2.487975); // As fitted by Cyclops
		sccs.makeDirty();
		System.err.println(sccs.getLogLikelihood());
		// Should be ~ -10.08828 according to Cyclops
	}
}
