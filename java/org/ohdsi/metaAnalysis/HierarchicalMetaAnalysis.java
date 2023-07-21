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
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.NormalDistribution;
import org.ohdsi.likelihood.CachedModelLikelihood;
import org.ohdsi.mcmc.Analysis;
import org.ohdsi.mcmc.Runner;
import org.ohdsi.simpleDesign.SimpleLinearModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class HierarchicalMetaAnalysis implements Analysis {

	private final Likelihood likelihood;
	private final Likelihood prior;
	private final Likelihood joint;
	private final List<Parameter> allParameters = new ArrayList<>();
	private final OperatorSchedule schedule;

	public HierarchicalMetaAnalysis(List<DataModel> allMetaAnalysisDataModels) {

		// TODO Pass as configuration
		double hierarchicalLocationHyperStdDev = 1.0;

		double gammaHyperScale = 1.0;
		double gammaHyperRate = 1.0;

		double exposureHyperLocation = 0.0;
		double exposureHyperStdDev = 2.0;

		double tauScale = 1.0;
		double tauRate = 1.0;

		AdaptationMode mode = AdaptationMode.ADAPTATION_ON;
		double operatorWeight = 1.0;

		double startingTau = 1.0;

		long seed = 666;

		MathUtils.setSeed(seed);

		// Build data likelihood, main-effect parameters and their operators
		List<Likelihood> allDataLikelihoods = new ArrayList<>();
		List<MCMCOperator> allOperators = new ArrayList<>();

		// Build data likelihood
		int metaAnalysisCount = 0;
		for (DataModel metaAnalysis : allMetaAnalysisDataModels) {

			allDataLikelihoods.add(
					new CachedModelLikelihood("likelihood" + (metaAnalysisCount + 1),
							metaAnalysis));  // TODO Do we need to cache? or just allLikelihoods.add(group.getLikelihood());

			Parameter beta = metaAnalysis.getCompoundParameter();
			allParameters.add(beta);

			allOperators.add(new RandomWalkOperator(beta, null, 0.75,
					RandomWalkOperator.BoundaryCondition.reflecting, operatorWeight * beta.getDimension(), mode)); // TODO Use HMC

			++metaAnalysisCount;
		}

		CompoundParameter allBetas = new CompoundParameter("all.beta");
		for (Parameter beta : allParameters) {
			allBetas.addParameter(beta);
		}

		this.likelihood = new CompoundLikelihood(allDataLikelihoods); // TODO Use multiple threads

		// Build meta-analysis model effects
		DesignMatrix designMatrix = new DesignMatrix("designMatrix",
				new Parameter[] { // TODO Rewrite for unbalanced dataModels
						// within-outcome bias
						makeIEffect("dm.outcome1",0, 3, 4),
						makeIEffect("dm.outcome2",1, 3, 4),
						makeIEffect("dm.outcome3",2, 3, 4),
						// within-data-source bias
						makeJEffect("dm.source1", 0, 3, 4),
						makeJEffect("dm.source2", 1, 3, 4),
						makeJEffect("dm.source3", 2, 3, 4),
						makeJEffect("dm.source4", 3, 3, 4),
						// de-biased exposure effect
						makeIEffect("dm.exposure", 2, 3, 4),
				}, false);

		Parameter outcomeEffect = randomize("outcome", 3, 0, 1); // TODO Scale for over-dispersed noise
		Parameter sourceEffect = randomize("source", 4, 0, 1);
		Parameter exposureEffect = randomize("exposure", 1, 0, 1);

		CompoundParameter allEffects = new CompoundParameter("allEffects",
				new Parameter[] {
						outcomeEffect,
						sourceEffect,
						exposureEffect,
				});

		if (designMatrix.getColumnDimension() != allEffects.getDimension()) {
			throw new RuntimeException("Invalid parameter dimensions");
		}

		// Build hierarchical priors and operators
		Parameter tau = new Parameter.Default("tau", startingTau, 0.0, Double.POSITIVE_INFINITY);

		SimpleLinearModel allEffectDistribution = new SimpleLinearModel("linearModel",
				allBetas, designMatrix, allEffects, tau);

		DistributionLikelihood tauPrior = new DistributionLikelihood(
				new GammaDistribution(tauScale, tauRate));
		tauPrior.addData(tau);

		MCMCOperator tauOperator = new ScaleOperator(tau, 0.75, mode, operatorWeight);

		List<Likelihood> allPriors = new ArrayList<>();

		allPriors.add(allEffectDistribution);
		allPriors.add(tauPrior);
		allParameters.add(tau);
		allOperators.add(tauOperator);

		HierarchicalNormalComponents outcomeComponents = makeHierarchicalNormalComponents(
				"outcome", outcomeEffect, new LocationHyperPrior(hierarchicalLocationHyperStdDev),
				new GammaOnPrecisionPrior(gammaHyperScale, gammaHyperRate), operatorWeight, mode);

		allParameters.add(outcomeEffect);
		allParameters.addAll(outcomeComponents.parameters);
		allOperators.addAll(outcomeComponents.operators);
		allPriors.addAll(outcomeComponents.likelihoods);

		HierarchicalNormalComponents sourceComponents = makeHierarchicalNormalComponents(
				"source", sourceEffect, new LocationHyperPrior(hierarchicalLocationHyperStdDev),
				new GammaOnPrecisionPrior(gammaHyperScale, gammaHyperScale), operatorWeight, mode);

		allParameters.add(sourceEffect);
		allParameters.addAll(sourceComponents.parameters);
		allOperators.addAll(sourceComponents.operators);
		allPriors.addAll(sourceComponents.likelihoods);

		// Build exposure prior and operator
		DistributionLikelihood exposureDistribution = new DistributionLikelihood(
				new NormalDistribution(exposureHyperLocation, exposureHyperStdDev));

		MCMCOperator exposureOperator = new RandomWalkOperator(exposureEffect, null, 0.75,
				RandomWalkOperator.BoundaryCondition.reflecting, operatorWeight, mode); // TODO Gibbs sample

		allParameters.add(exposureEffect);
		allOperators.add(exposureOperator);
		allPriors.add(exposureDistribution);

		this.prior = new CompoundLikelihood(allPriors); // TODO Use multiple threads
		this.joint = new CompoundLikelihood(Arrays.asList(likelihood, prior)); // TODO Use multiple threads
		this.joint.setId("joint");

		this.schedule = new SimpleOperatorSchedule(1000, 0.0);
		this.schedule.addOperators(allOperators);
	}

	@Override
	 public List<Loggable> getLoggerColumns() {

		List<Loggable> columns = new ArrayList<>();
		columns.add(likelihood);
		columns.add(prior);
		columns.addAll(allParameters);

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

	public static Parameter makeIEffect(String name, int index, int I, int J) {
		double[] effect = new double[I * J];
		for (int j = 0; j < J; ++j) {
			effect[index * J + j] = 1.0;
		}
		return new Parameter.Default(name, effect);
	}

	public static Parameter makeJEffect(String name, int index, int I, int J) {
		double[] effect = new double[I * J];
		for (int i = 0; i < I; ++i) {
			effect[i * J + index] = 1.0;
		}
		return new Parameter.Default(name, effect);
	}

	public static Parameter randomize(String name, int dim, double center, double scale) {
		double[] effect = new double[dim];
		for (int i = 0; i < dim; ++i) {
			effect[i] = center + scale * MathUtils.nextGaussian();
		}
		Parameter parameter = new Parameter.Default(name, effect);
		parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, dim));
		return parameter;
	}

	static class HierarchicalNormalComponents {

		final List<Parameter> parameters;
		final List<MCMCOperator> operators;
		final List<Likelihood> likelihoods;

		HierarchicalNormalComponents(List<Parameter> parameters,
									 List<MCMCOperator> operators,
									 List<Likelihood> likelihoods) {
			this.parameters = parameters;
			this.operators = operators;
			this.likelihoods = likelihoods;
		}
	}

	public static HierarchicalNormalComponents makeHierarchicalNormalComponents(String name,
																				Parameter effects,
																				LocationHyperPrior locationHyperPrior,
																				ScalePrior scalePrior, double weight,
																				AdaptationMode mode) {

		Parameter mean = randomize(name + ".mean", 1, 0, 1);
		Parameter scale = scalePrior.getParameter();
		scale.setId(name + ".precision");

		DistributionLikelihood distribution = new DistributionLikelihood(
				new NormalDistributionModel(mean, scale, scalePrior.isPrecision()));
		distribution.addData(effects);

		DistributionLikelihood meanHyperDistribution = new DistributionLikelihood(
				new NormalDistribution(0.0, locationHyperPrior.getSd()));
		meanHyperDistribution.addData(mean);

		List<Parameter> parameters = new ArrayList<>();
		parameters.add(mean);
		parameters.add(scale);

		List<MCMCOperator> operators = new ArrayList<>();
		operators.add(new RandomWalkOperator(effects, null, 0.75,
				RandomWalkOperator.BoundaryCondition.reflecting, weight * effects.getDimension(), mode)); // TODO Gibbs sample
		operators.add(new RandomWalkOperator(mean, null, 0.75,
				RandomWalkOperator.BoundaryCondition.reflecting, weight, mode)); // TODO Gibbs sample
		operators.add(scalePrior.getOperator(distribution, weight, mode));

		List<Likelihood> likelihood = new ArrayList<>();
		likelihood.add(distribution);
		likelihood.add(meanHyperDistribution);
		likelihood.add(scalePrior.getPrior());

		return new HierarchicalNormalComponents(parameters, operators, likelihood);
	}

	public static void main(String[] args) {

		int chainLength = 1100000;
		int burnIn = 100000;
		int subSampleFrequency = 1000;


		List<DataModel> allDataModels = new ArrayList<>();
		allDataModels.add(new ExtendingEmpiricalDataModel("ForDavid/grids_example_1.csv"));
		allDataModels.add(new ExtendingEmpiricalDataModel("ForDavid/grids_example_2.csv"));
		allDataModels.add(new ExtendingEmpiricalDataModel("ForDavid/grids_example_3.csv"));

		HierarchicalMetaAnalysis analysis = new HierarchicalMetaAnalysis(allDataModels);

		Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, 666);

		runner.run();

		runner.processSamples();
	}
}
