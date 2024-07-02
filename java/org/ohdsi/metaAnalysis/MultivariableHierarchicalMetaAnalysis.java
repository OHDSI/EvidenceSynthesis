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
package org.ohdsi.metaAnalysis;

import dr.evomodel.operators.PrecisionMatrixGibbsOperator;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.math.distributions.*;
import dr.math.matrixAlgebra.IllegalDimension;
import org.ohdsi.likelihood.ConditionalPoissonLikelihood;
import org.ohdsi.likelihood.MultivariableCoxPartialLikelihood;
import org.ohdsi.mcmc.Analysis;
import org.ohdsi.mcmc.Runner;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.ohdsi.likelihood.MultivariableCoxPartialLikelihood.exampleBladder;

public class MultivariableHierarchicalMetaAnalysis implements Analysis {

	private final Likelihood likelihood;
	private final Likelihood prior;
	private final Likelihood joint;
	private final List<Parameter> parameters;
	private final OperatorSchedule schedule;

	public MultivariableHierarchicalMetaAnalysis(List<ConditionalPoissonLikelihood> likelihoods,
												 HierarchicalMetaAnalysisConfiguration cg) {

		MathUtils.setSeed(cg.seed);

		// Build data likelihood, main-effect parameters and their operators
		List<Likelihood> allDataLikelihoods = new ArrayList<>();
		List<MCMCOperator> allOperators = new ArrayList<>();
		List<Parameter> allParameters = new ArrayList<>();

		final int analysisDim = likelihoods.get(0).getParameter().getDimension();

		// Build data likelihood
		for (ConditionalPoissonLikelihood singleAnalysis : likelihoods) {

			if (analysisDim != singleAnalysis.getParameter().getDimension()) {
				throw new IllegalArgumentException("Mismatched regression dimensions");
			}

			allDataLikelihoods.add(singleAnalysis);
			Parameter beta = singleAnalysis.getParameter();
			allParameters.add(beta);

			allOperators.add(new RandomWalkOperator(beta, null, 0.1, // TODO HMC will be way faster!!!
					RandomWalkOperator.BoundaryCondition.reflecting, cg.operatorWeight * beta.getDimension(), cg.mode));
		}
		// End of data likelihood

		// Build hierarchical priors and operators
		Parameter mu = new Parameter.Default("mean", analysisDim, cg.startingMu);
		mu.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, analysisDim));

//		Parameter tau = new Parameter.Default("tau", analysisDim, cg.startingTau);
//		tau.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, analysisDim));
		MatrixParameter tau2 = diagonalMatrixParameter("tau", analysisDim, cg.startingTau);

		MultivariateDistributionLikelihood hierarchy = new MultivariateDistributionLikelihood(
			new MultivariateNormalDistributionModel(mu,
//					new DiagonalMatrix(tau)
					tau2
			));

		for (Parameter beta : allParameters) {
			hierarchy.addData(beta);
		}

		double[] muPriorMean = new double[analysisDim];
		Arrays.fill(muPriorMean, cg.muMean);
		double muPriorPrecision = 1 / (cg.muSd * cg.muSd);

		MultivariateDistributionLikelihood meanPrior = new MultivariateDistributionLikelihood(
			new MultivariateNormalDistribution(muPriorMean,muPriorPrecision));
		meanPrior.addData(mu);

//		DistributionLikelihood tauPrior = new DistributionLikelihood(
//				new GammaDistribution(cg.tauShape, cg.tauScale));
//		tauPrior.addData(tau);

		MultivariateDistributionLikelihood tau2Prior = new MultivariateDistributionLikelihood(
				new WishartDistribution(cg.tauDf, diagonalScaleMatrix(analysisDim, cg.tauScale))
		);
		tau2Prior.addData(tau2);

//		List<Likelihood> allPriors = Arrays.asList(hierarchy, meanPrior, tauPrior);
		List<Likelihood> allPriors = Arrays.asList(hierarchy, meanPrior, tau2Prior);

		allParameters.add(mu);
//		allParameters.add(tau);
		allParameters.add(tau2);


		MCMCOperator meanOperator = null;
		try {
			meanOperator = new MultivariateNormalGibbsOperator(hierarchy, meanPrior, 1.0);
		} catch (IllegalDimension e) {
			e.printStackTrace();
		}

//		MCMCOperator tauOperator = new ScaleOperator(tau, 0.75, cg.mode, 1.0);
		MCMCOperator tau2Operator = new PrecisionMatrixGibbsOperator(hierarchy,
				(WishartStatistics) tau2Prior.getDistribution(), 1.0);

		allOperators.add(meanOperator);
//		allOperators.add(tauOperator);
		allOperators.add(tau2Operator);


		// Finalize
		this.prior = new CompoundLikelihood(allPriors);
		this.likelihood = new CompoundLikelihood(cg.threads, allDataLikelihoods);
		this.joint = new CompoundLikelihood(Arrays.asList(likelihood, prior));
		this.joint.setId("joint");

		this.parameters = allParameters;
		this.schedule = new SimpleOperatorSchedule(1000, 0.0);
		this.schedule.addOperators(allOperators);
	}

	@Override
	 public List<Loggable> getLoggerColumns() {

		List<Loggable> columns = new ArrayList<>();
		columns.add(likelihood);
		columns.add(prior);
		columns.addAll(parameters);

		return columns;
	}

	@Override
	public Likelihood getJoint() {
		return joint;
	}

	@Override
	public OperatorSchedule getSchedule() {
		return schedule;
	}

	public static class HierarchicalMetaAnalysisConfiguration {

		// make all fields public to allow rJava interface
		// alternatively, can also do: (a) turn into constructor (but no default allowed); (b). write a setter for each field (too much trouble)

		// gamma prior for hierarchy precision
		public double tauShape = 1;
		public double tauScale = 1;
		public double tauDf = 3;

		// normal prior for hierarchy mean
		public double muMean = 0;
		public double muSd = 2;

		public double startingMu = 0;
		public double startingTau = 1;

		AdaptationMode mode = AdaptationMode.ADAPTATION_ON;
		public double operatorWeight = 10.0;
		public long seed = 666;

		public int threads = 1;
	}


	public static void main(String[] args) {

		int chainLength = 110000;
		int burnIn = 10000;
		int subSampleFrequency = 10;

		HierarchicalMetaAnalysisConfiguration cg = new HierarchicalMetaAnalysisConfiguration();

		List<ConditionalPoissonLikelihood> likelihoods = Arrays.asList(
				new MultivariableCoxPartialLikelihood(
						new Parameter.Default(new double[] { -0.4608773, -0.1012988 }),
						exampleBladder()),
				new MultivariableCoxPartialLikelihood(
						new Parameter.Default(new double[] { -0.4608773, -0.1012988 }),
						exampleBladder()),
				new MultivariableCoxPartialLikelihood(
						new Parameter.Default(new double[] { -0.4608773, -0.1012988 }),
						exampleBladder()),
				new MultivariableCoxPartialLikelihood(
						new Parameter.Default(new double[] { -0.4608773, -0.1012988 }),
						exampleBladder())
		);

		MultivariableHierarchicalMetaAnalysis analysis = new MultivariableHierarchicalMetaAnalysis(likelihoods, cg);

		Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, cg.seed);

		runner.run();

		runner.processSamples();
	}

	public static double[][] diagonalScaleMatrix(int dim, double diagonal) {
		double[][] scale = new double[dim][];
		for (int i = 0; i < dim; ++i) {
			scale[i] = new double[dim];
			scale[i][i] = diagonal;
		}
		return scale;
	}

	public static MatrixParameter diagonalMatrixParameter(String name, int dim, double diagonal) {
		Parameter[] columns = new Parameter[dim];
		for (int i = 0; i < dim; ++i) {
			columns[i] = new Parameter.Default(dim, 0.0);
			columns[i].setParameterValue(i, diagonal);
		}
		return new MatrixParameter(name, columns);
	}
}
