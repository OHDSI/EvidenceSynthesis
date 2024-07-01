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

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.NormalDistribution;
import org.ohdsi.likelihood.ConditionalPoissonLikelihood;
import org.ohdsi.mcmc.Analysis;
import org.ohdsi.mcmc.Runner;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MultivariableHierarchicalMetaAnalysis implements Analysis {

	private final Likelihood likelihood;
	private final Likelihood prior;
	private final Likelihood joint;
	private final List<Parameter> parameters;
	private final OperatorSchedule schedule;

	public MultivariableHierarchicalMetaAnalysis(List<ConditionalPoissonLikelihood> likelihoods,
												 HierarchicalMetaAnalysisConfiguration cg) {

		MathUtils.setSeed(666);

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

			allOperators.add(new RandomWalkOperator(beta, null, 0.75, // TODO HMC will be way faster!!!
					RandomWalkOperator.BoundaryCondition.reflecting, cg.operatorWeight * beta.getDimension(), cg.mode)); // TODO Use HMC
		}
		// End of data likelihood

		// Build hierarchical priors and operators
		Parameter mu = new Parameter.Default("mean", analysisDim, cg.startingMu);
		Parameter tau = new Parameter.Default("tau", analysisDim, cg.startingTau);
		tau.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, analysisDim));

		MultivariateDistributionLikelihood hierarchy = new MultivariateDistributionLikelihood(
			new MultivariateNormalDistributionModel(mu, new DiagonalMatrix(tau)));

		for (Parameter beta : allParameters) {
			hierarchy.addData(beta);
		}

		DistributionLikelihood meanPrior = new DistributionLikelihood(
				new NormalDistribution(cg.muMean, cg.muSd));

		DistributionLikelihood tauPrior = new DistributionLikelihood(
				new GammaDistribution(cg.tauShape, cg.tauScale));

		List<Likelihood> allPriors = Arrays.asList(hierarchy, meanPrior, tauPrior);
		allParameters.add(mu);
		allParameters.add(tau);

		MCMCOperator meanOperator = new RandomWalkOperator(mu, null, 0.75,
				RandomWalkOperator.BoundaryCondition.reflecting, cg.operatorWeight, cg.mode);
		MCMCOperator tauOperator = new ScaleOperator(tau, 0.75, cg.mode, cg.operatorWeight);
		allOperators.add(meanOperator);
		allOperators.add(tauOperator);

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

		// gamma prior for std of the random error
		public double tauShape = 1.0;
		public double tauScale = 1.0;

		public double startingMu = 0.0;
		public double startingTau = 1.0;

		AdaptationMode mode = AdaptationMode.ADAPTATION_ON;
		public double operatorWeight = 1.0;

		public long seed = 666;

		public int threads = 1;

		public double muMean = 0;
		public double muSd = 10;
	}


	public static void main(String[] args) {

		int chainLength = 1100000;
		int burnIn = 100000;
		int subSampleFrequency = 1000;

		HierarchicalMetaAnalysisConfiguration cg = new HierarchicalMetaAnalysisConfiguration();

		MultivariableHierarchicalMetaAnalysis analysis = new MultivariableHierarchicalMetaAnalysis(null, cg);

		Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, cg.seed);

		runner.run();

		runner.processSamples();
	}
}
