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

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.hmc.CompoundGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.JointGradient;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.inference.operators.*;
import dr.inference.operators.hmc.*;
import dr.inference.regression.CyclopsRegressionModelGradient;
import dr.math.MathUtils;
import org.ohdsi.likelihood.MultivariableCoxPartialLikelihood;
import org.ohdsi.mcmc.Analysis;
import org.ohdsi.mcmc.Runner;

import java.sql.Array;
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

	public MultivariableHierarchicalMetaAnalysis(List<DataModel> dataModels,
												 MultivariatePrior multivariatePrior,
												 HierarchicalMetaAnalysisConfiguration cg) {

		MathUtils.setSeed(cg.seed);

		// Build data likelihood, main-effect parameters and their operators
		List<Likelihood> allDataLikelihoods = new ArrayList<>();
		List<MCMCOperator> allOperators = new ArrayList<>();
		List<Parameter> allParameters = new ArrayList<>();

		final int analysisDim = dataModels.get(0).getCompoundParameter().getDimension();

		List<GradientWrtParameterProvider> likelihoodDerivativeList = new ArrayList<>();
		List<GradientWrtParameterProvider> priorDerivativeList = new ArrayList<>();
		final MultivariateDistributionLikelihood mdl = (MultivariateDistributionLikelihood) multivariatePrior.getLikelihood(0);
		final GradientProvider provider = (GradientProvider) mdl.getDistribution();
//		List<Parameter> allBetas = new ArrayList<>();

		// Build data likelihood
		for (DataModel singleAnalysis : dataModels) {

			if (analysisDim != singleAnalysis.getCompoundParameter().getDimension()) {
				throw new IllegalArgumentException("Mismatched regression dimensions");
			}

			allDataLikelihoods.add(singleAnalysis.getLikelihood());
			Parameter beta = singleAnalysis.getCompoundParameter();
//			allBetas.add(beta);
//			allParameters.add(beta);

//			allOperators.add(new RandomWalkOperator(beta, null, 0.1, // TODO HMC will be way faster!!!
//					RandomWalkOperator.BoundaryCondition.reflecting, cg.operatorWeight * beta.getDimension(), cg.mode));
			likelihoodDerivativeList.add((GradientWrtParameterProvider) singleAnalysis.getLikelihood());
			priorDerivativeList.add(new GradientWrtParameterProvider.ParameterWrapper(provider, beta, mdl));
//			allOperators.add(new HamiltonianMonteCarloOperator(AdaptationMode.ADAPTATION_OFF, 0.1, ))

		}

		final CompoundGradient priorGradient = new CompoundGradient(priorDerivativeList);
		final CompoundGradient likelihoodGradient = new CompoundGradient(likelihoodDerivativeList);

		allParameters.add(priorGradient.getParameter());

		List<GradientWrtParameterProvider> jointGradientList = Arrays.asList(priorGradient, likelihoodGradient);
		JointGradient jointGradient = new JointGradient(jointGradientList);
		jointGradient.getGradientLogDensity();

		final double stepSize = 1.8;
		final int nSteps = 10;
		final double randomStepFraction = 0;
		final MassPreconditioningOptions preconditioningOptions =
				new MassPreconditioningOptions.Default(10, 0, 0, 0, false, new Parameter.Default(1E-2), new Parameter.Default(1E2));
		final MassPreconditioner.Type preconditionerType = MassPreconditioner.Type.DIAGONAL;
		final MassPreconditioner preconditioner = preconditionerType.factory(jointGradient, null, preconditioningOptions);

		HamiltonianMonteCarloOperator.Options runtimeOptions = new HamiltonianMonteCarloOperator.Options(
				stepSize, nSteps, randomStepFraction,
				preconditioningOptions,
				0, 0,
				10, 0.1,
				0.8,
				HamiltonianMonteCarloOperator.InstabilityHandler.factory("reject"));

		allOperators.add(new HamiltonianMonteCarloOperator(AdaptationMode.ADAPTATION_ON, 0.5, jointGradient, jointGradient.getParameter(), null, null, runtimeOptions, preconditioner));

		// End of data likelihood

		// Build hierarchical priors and operators
//		Parameter mu = new Parameter.Default("mean", analysisDim, cg.startingMu);
//		mu.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, analysisDim));
//
////		Parameter tau = new Parameter.Default("tau", analysisDim, cg.startingTau);
////		tau.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, analysisDim));
//		MatrixParameter tau2 = diagonalMatrixParameter("tau", analysisDim, cg.startingTau);
//
//		MultivariateDistributionLikelihood hierarchy = new MultivariateDistributionLikelihood(
//			new MultivariateNormalDistributionModel(mu,
////					new DiagonalMatrix(tau)
//					tau2
//			));
//
//		for (Parameter beta : allParameters) {
//			hierarchy.addData(beta);
//		}
//
//		double[] muPriorMean = new double[analysisDim];
//		Arrays.fill(muPriorMean, cg.muMean);
//		double muPriorPrecision = 1 / (cg.muSd * cg.muSd);
//
//		MultivariateDistributionLikelihood meanPrior = new MultivariateDistributionLikelihood(
//			new MultivariateNormalDistribution(muPriorMean,muPriorPrecision));
//		meanPrior.addData(mu);
//
////		DistributionLikelihood tauPrior = new DistributionLikelihood(
////				new GammaDistribution(cg.tauShape, cg.tauScale));
////		tauPrior.addData(tau);
//
//		MultivariateDistributionLikelihood tau2Prior = new MultivariateDistributionLikelihood(
//				new WishartDistribution(cg.tauDf, diagonalScaleMatrix(analysisDim, cg.tauScale))
//		);
//		tau2Prior.addData(tau2);
//
////		List<Likelihood> allPriors = Arrays.asList(hierarchy, meanPrior, tauPrior);
//		List<Likelihood> allPriors = Arrays.asList(hierarchy, meanPrior, tau2Prior);
//
//		allParameters.add(mu);
////		allParameters.add(tau);
//		allParameters.add(tau2);
//
//
//		MCMCOperator meanOperator = null;
//		try {
//			meanOperator = new MultivariateNormalGibbsOperator(hierarchy, meanPrior, 1.0);
//		} catch (IllegalDimension e) {
//			e.printStackTrace();
//		}
//
////		MCMCOperator tauOperator = new ScaleOperator(tau, 0.75, cg.mode, 1.0);
//		MCMCOperator tau2Operator = new PrecisionMatrixGibbsOperator(hierarchy,
//				(WishartStatistics) tau2Prior.getDistribution(), 1.0);
//
//		allOperators.add(meanOperator);
////		allOperators.add(tauOperator);
//		allOperators.add(tau2Operator);

		allParameters.addAll(multivariatePrior.getParameters());
		allOperators.addAll(multivariatePrior.getOperators(cg.operatorWeight, cg.mode));

		// Finalize
//		this.prior = new CompoundLikelihood(allPriors);
		this.prior = multivariatePrior.getPrior();
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

		List<DataModel> likelihoods = Arrays.asList(
				new MultivariableCoxPartialLikelihood(
						new Parameter.Default(new double[] { -0.4608773, -0.1012988 }),
						exampleBladder()),
				new MultivariableCoxPartialLikelihood(
						new Parameter.Default(new double[] { -0.5608773, -0.2012988 }),
						exampleBladder()),
				new MultivariableCoxPartialLikelihood(
						new Parameter.Default(new double[] { -0.6608773, -0.3012988 }),
						exampleBladder()),
				new MultivariableCoxPartialLikelihood(
						new Parameter.Default(new double[] { -0.7608773, -0.4012988 }),
						exampleBladder())
		);

		MultivariableHierarchicalMetaAnalysis analysis = new MultivariableHierarchicalMetaAnalysis(likelihoods,
				new MultivariatePrior.MultivariateNormal(likelihoods, cg), cg);

		Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, cg.seed);

		runner.run();

		runner.processSamples();

		System.exit(0); // Close threads
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
