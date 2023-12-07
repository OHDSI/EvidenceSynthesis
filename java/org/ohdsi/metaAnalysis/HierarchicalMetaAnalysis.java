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
import dr.inference.hmc.CompoundDerivative;
import dr.inference.hmc.CompoundGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.JointGradient;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.inference.operators.*;
import dr.inference.operators.hmc.*;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.NormalDistribution;
import org.ohdsi.hmc.HmcOptions;
import org.ohdsi.likelihood.CachedModelLikelihood;
import org.ohdsi.mcmc.Analysis;
import org.ohdsi.mcmc.Runner;
import org.ohdsi.simpleDesign.SimpleLinearModel;
import org.ohdsi.simpleDesign.SimpleLinearModelGradientWrtArgument;
import org.ohdsi.simpleDesign.SimpleLinearModelGradientWrtEffects;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class HierarchicalMetaAnalysis implements Analysis {

	private final Likelihood likelihood;
	private final Likelihood prior;
	private final Likelihood joint;
	private final List<Parameter> parameters;
	private final OperatorSchedule schedule;

	private ScalePrior makeHierarchicalScalePrior(double hyper1, double hyper2) {
		return new GammaOnPrecisionPrior(hyper1, hyper2);
	}

	public HierarchicalMetaAnalysis(List<DataModel> allMetaAnalysisDataModels,
									HierarchicalMetaAnalysisConfiguration cg) {

		MathUtils.setSeed(cg.seed);

		// Build data likelihood, main-effect parameters and their operators
		List<Likelihood> allDataLikelihoods = new ArrayList<>();
		List<MCMCOperator> allOperators = new ArrayList<>();
		List<Parameter> allParameters = new ArrayList<>();

		// Build data likelihood
		int metaAnalysisCount = 0;
		for (DataModel metaAnalysis : allMetaAnalysisDataModels) {

			int maLabel = (metaAnalysisCount + 1);
			allDataLikelihoods.add(
					new CachedModelLikelihood("likelihood" + maLabel,
							metaAnalysis));  // TODO Do we need to cache? or just allLikelihoods.add(group.getLikelihood());

			Parameter beta = metaAnalysis.getCompoundParameter();
			allParameters.add(beta);

			allOperators.add(new RandomWalkOperator(beta, null, 0.75,
					RandomWalkOperator.BoundaryCondition.reflecting, cg.operatorWeight * beta.getDimension(), cg.mode)); // TODO Use HMC

			++metaAnalysisCount;
		}

		CompoundParameter allBetas = new CompoundParameter("all.beta");
		for (Parameter beta : allParameters) {
			allBetas.addParameter(beta);
		}
		// End of data likelihood

		// Build hierarchical priors and operators
		DesignMatrix designMatrix = new DesignMatrix("designMatrix", false);
		CompoundParameter allEffects = new CompoundParameter("allEffects");
		List<Likelihood> allPriors = new ArrayList<>();
		List<GradientWrtParameterProvider> allEffectsGradient = new ArrayList<>();

		Parameter tau = new Parameter.Default("tau", cg.startingTau, 0.0, Double.POSITIVE_INFINITY);
		DistributionLikelihood tauPrior = new DistributionLikelihood(
				new GammaDistribution(cg.tauShape, cg.tauScale));
		tauPrior.addData(tau);
		MCMCOperator tauOperator = new ScaleOperator(tau, 0.75, cg.mode, cg.operatorWeight);
		allPriors.add(tauPrior);
		allParameters.add(tau);
		allOperators.add(tauOperator);

		int primaryCount = addPrimaryDesign(designMatrix, allMetaAnalysisDataModels, cg.primaryEffectName, cg.separateEffectPrior);
		Parameter primaryEffect = randomize(cg.primaryEffectName, primaryCount, 0, 1); // TODO Scale for over-dispersed noise
		allEffects.addParameter(primaryEffect);

		HierarchicalNormalComponents primaryComponents = makeHierarchicalNormalComponents(
				cg.primaryEffectName, primaryEffect, cg.hierarchicalLocationPrimaryHyperStdDev,
				makeHierarchicalScalePrior(cg.gammaHyperPrimaryShape, cg.gammaHyperPrimaryScale), cg.operatorWeight, cg.mode);

		allParameters.add(primaryEffect);
		allParameters.addAll(primaryComponents.parameters);
		allOperators.addAll(primaryComponents.operators);
		allPriors.addAll(primaryComponents.likelihoods);
		allEffectsGradient.addAll(primaryComponents.gradients);

		if (cg.includeSecondary) {

			int secondaryCount = addSecondaryDesign(designMatrix, allMetaAnalysisDataModels, cg.secondaryEffectName);
			Parameter secondaryEffect = randomize(cg.secondaryEffectName, secondaryCount, 0, 1);
			allEffects.addParameter(secondaryEffect);

			HierarchicalNormalComponents secondaryComponents = makeHierarchicalNormalComponents(
					cg.secondaryEffectName, secondaryEffect, cg.hierarchicalLocationSecondaryHyperStdDev,
					makeHierarchicalScalePrior(cg.gammaHyperSecondaryShape, cg.gammaHyperSecondaryScale), cg.operatorWeight, cg.mode);

			allParameters.add(secondaryEffect);
			allParameters.addAll(secondaryComponents.parameters);
			allOperators.addAll(secondaryComponents.operators);
			allPriors.addAll(secondaryComponents.likelihoods);
			allEffectsGradient.addAll(secondaryComponents.gradients);
		}

		if (cg.includeExposure) {
			int effectCount = addEffectDesign(designMatrix, allMetaAnalysisDataModels, cg.exposureEffectName);
			Parameter exposureEffect = randomize(cg.exposureEffectName, effectCount, 0, 1);
			allEffects.addParameter(exposureEffect);

			DistributionLikelihood exposureDistribution = new DistributionLikelihood(
					new NormalDistribution(cg.exposureHyperLocation, cg.exposureHyperStdDev));
			exposureDistribution.addData(exposureEffect);

			MCMCOperator exposureOperator = new RandomWalkOperator(exposureEffect, null, 0.75,
					RandomWalkOperator.BoundaryCondition.reflecting, cg.operatorWeight, cg.mode); // TODO Gibbs sample

			allParameters.add(exposureEffect);
			allOperators.add(exposureOperator);
			allPriors.add(exposureDistribution);
			if (exposureDistribution.getDistribution() instanceof GradientProvider) {
				allEffectsGradient.add(new GradientWrtParameterProvider.ParameterWrapper(
						(GradientProvider) exposureDistribution.getDistribution(), exposureEffect, exposureDistribution));
			}

			//System.err.println("added exposure effect!");
		}

		if (designMatrix.getColumnDimension() != allEffects.getDimension()) {
			throw new RuntimeException("Invalid parameter dimensions");
		}

		SimpleLinearModel allEffectDistribution = new SimpleLinearModel("linearModel",
				allBetas, designMatrix, allEffects, tau);
		allPriors.add(allEffectDistribution);

		// Finalize
		this.prior = new CompoundLikelihood(allPriors); // TODO Use multiple threads?
		this.likelihood = new CompoundLikelihood(allDataLikelihoods); // TODO Use multiple threads?
		this.joint = new CompoundLikelihood(Arrays.asList(likelihood, prior));
		this.joint.setId("joint");

		this.parameters = allParameters;

//		GradientWrtParameterProvider subGradientDataModel1 = makeDataModelCompoundGradient(allMetaAnalysisDataModels);
//		GradientWrtParameterProvider subGradientDataModel2 = new SimpleLinearModelGradientWrtArgument(allEffectDistribution);
//		JointGradient gradientDataModel = new JointGradient(List.of(subGradientDataModel1, subGradientDataModel2));
//
//		System.err.println(gradientDataModel.getReport());
//
//		GradientWrtParameterProvider subGradientEffects1 = new SimpleLinearModelGradientWrtEffects(allEffectDistribution);
//		CompoundGradient subGradientEffects2 = new CompoundDerivative(allEffectsGradient);
//		JointGradient gradientEffects = new JointGradient(List.of(subGradientEffects1, subGradientEffects2));
//
//		System.err.println(gradientEffects.getReport());
//
//		CompoundGradient totalGradient = new CompoundDerivative(List.of(gradientDataModel, gradientEffects));
//
//		System.err.println(totalGradient.getReport());
//
//		// Use HMC
//		HmcOptions options = new HmcOptions(totalGradient);
//		HamiltonianMonteCarloOperator hmc = new HamiltonianMonteCarloOperator(
//				AdaptationMode.ADAPTATION_ON, 1.0,
//				totalGradient,
//				totalGradient.getParameter(), null, null,
//				options.getHmcOptions(), options.getPreconditioner());
//
//		// Remove old MH operators
//		allOperators.remove(0);
//		allOperators.remove(0);
//		allOperators.remove(0);
//		allOperators.remove(1);
//		allOperators.remove(3);
//		allOperators.remove(5);
//
//		allOperators.add(hmc);

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

	public static Parameter randomize(String name, int dim, double center, double scale) {
		double[] effect = new double[dim];
		for (int i = 0; i < dim; ++i) {
			effect[i] = center + scale * MathUtils.nextGaussian();
		}
		Parameter parameter = new Parameter.Default(name, effect);
		parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, dim));
		return parameter;
	}

	public static class HierarchicalMetaAnalysisConfiguration {

		// make all fields public to allow rJava interface
		// alternatively, can also do: (a) turn into constructor (but no default allowed); (b). write a setter for each field (too much trouble)

		//prior standard deviation for primary & secondary effect mean
		public double hierarchicalLocationPrimaryHyperStdDev = 1.0;
		public double hierarchicalLocationSecondaryHyperStdDev = 1.0;

		// gamma prior for primary & secondary effect precision (normal dist)
		public double gammaHyperPrimaryShape = 1.0;
		public double gammaHyperPrimaryScale = 1.0;

		public double gammaHyperSecondaryShape = 1.0;
		public double gammaHyperSecondaryScale = 1.0;

		// prior mean and std for exposure effect (global effect for outcome of interest)
		public double exposureHyperLocation = 0.0;
		public double exposureHyperStdDev = 2.0;

		// gamma prior for std of the random error
		public double tauShape = 1.0;
		public double tauScale = 1.0;

		public double startingTau = 1.0;

		AdaptationMode mode = AdaptationMode.ADAPTATION_ON;
		public double operatorWeight = 1.0;

		public long seed = 666;

		public String primaryEffectName = "outcome";
		public String secondaryEffectName = "source";
		public String exposureEffectName = "exposure";

		// include the secondary effects? (e.g., data source effects)
		public boolean includeSecondary = true;

		// include the exposure effect for main outcome of interest? (in addition to negative controls)
		public boolean includeExposure = true;

		// if using a separate prior on the main effect (i.e., set prior on biased effect instead of true effect)?
		public boolean separateEffectPrior = false;
	}

	static class HierarchicalNormalComponents {

		final List<Parameter> parameters;
		final List<MCMCOperator> operators;
		final List<Likelihood> likelihoods;
		final List<GradientWrtParameterProvider> gradients;

		HierarchicalNormalComponents(List<Parameter> parameters,
									 List<MCMCOperator> operators,
									 List<Likelihood> likelihoods,
									 List<GradientWrtParameterProvider> gradients) {
			this.parameters = parameters;
			this.operators = operators;
			this.likelihoods = likelihoods;
			this.gradients = gradients;
		}
	}

	public static HierarchicalNormalComponents makeHierarchicalNormalComponents(String name,
																				Parameter effects,
																				double locationHyperStdDev,
																				ScalePrior scalePrior, double weight,
																				AdaptationMode mode) {

		Parameter mean = randomize(name + ".mean", 1, 0, 1);
		Parameter scale = scalePrior.getParameter();
		scale.setId(name + ".scale");

		DistributionLikelihood distribution = new DistributionLikelihood(
				new NormalDistributionModel(mean, scale, scalePrior.isPrecision()));
		distribution.addData(effects);

		DistributionLikelihood meanHyperDistribution = new DistributionLikelihood(
				new NormalDistribution(0.0, locationHyperStdDev));
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

		List<GradientWrtParameterProvider> gradients = new ArrayList<>();
		if (distribution.getDistribution() instanceof GradientProvider) {
			gradients.add(new GradientWrtParameterProvider.ParameterWrapper(
					(GradientProvider) distribution.getDistribution(),
					effects, distribution));
		}

		return new HierarchicalNormalComponents(parameters, operators, likelihood, gradients);
	}

	private int addPrimaryDesign(DesignMatrix designMatrix,
								 List<DataModel> dataModels,
								 String effectName,
								 boolean separateEffectPrior) {
		int totalLength = getTotalNumberOfProfiles(dataModels);
		int effectOffset = totalLength;
		if (separateEffectPrior) {
			effectOffset -= dataModels.get(dataModels.size() - 1).getCompoundParameter().getDimension();
		}

		int offset = 0;
		int label = 0;
		for (DataModel dataModel : dataModels) {
			int length = dataModel.getCompoundParameter().getDimension();
			double[] effect = new double[totalLength];
			if (offset < effectOffset) {
				for (int i = 0; i < length; ++i) {
					effect[offset + i] = 1.0;
				}
			}
			designMatrix.addParameter(
					new Parameter.Default("dm." + effectName + (label + 1), effect));

			++label;
			offset += length;

		}
		return label;
	}

	private int addSecondaryDesign(DesignMatrix designMatrix,
								   List<DataModel> dataModels,
								   String effectName) {
		int totalLength = getTotalNumberOfProfiles(dataModels);
		int maxIdentifier = getMaxIdentifier(dataModels);

		for (int id = 0; id < maxIdentifier; ++id) {
			double[] effect = new double[totalLength];

			int offset = 0;
			for (DataModel dataModel : dataModels) {
				int length = dataModel.getCompoundParameter().getDimension();
				int whichIndex = findIdentifier(dataModel, id + 1);
				if (whichIndex >= 0) {
					effect[offset + whichIndex] = 1.0;
					// matching secondary (e.g., data source) effects
				}
				offset += length;
			}

			designMatrix.addParameter(
					new Parameter.Default("dm." + effectName + (id + 1), effect));
		}
		return maxIdentifier;
	}

	public int addEffectDesign(DesignMatrix designMatrix,
								List<DataModel> dataModels,
								String effectName) {
		int totalLength = getTotalNumberOfProfiles(dataModels);
		int offset = totalLength - dataModels.get(dataModels.size() - 1).getCompoundParameter().getDimension();

		double[] effect = new double[totalLength];
		for (int i = offset; i < totalLength; ++i) {
			effect[i] = 1.0;
		}

		designMatrix.addParameter(
				new Parameter.Default("dm." + effectName, effect));
		return 1;
	}

	private int getTotalNumberOfProfiles(List<DataModel> dataModels) {
		int length = 0;
		for (DataModel dataModel : dataModels) {
			length += dataModel.getCompoundParameter().getDimension();
		}
		return length;
	}

	private int getMaxIdentifier(List<DataModel> dataModels) {
		int max = maxOfList(dataModels.get(0).getIdentifiers());
		for (int i = 1; i < dataModels.size(); ++i) {
			max = Math.max(max, maxOfList(dataModels.get(i).getIdentifiers()));
		}
		return max;
	}

	private int maxOfList(List<Integer> integers) {
		int max = integers.get(0);
		for (int i = 1; i < integers.size(); ++i) {
			max = Math.max(max, integers.get(i));
		}
		return max;
	}

	private int findIdentifier(DataModel dataModel, int id) {
		List<Integer> identifiers = dataModel.getIdentifiers();
		return identifiers.indexOf(id);
	}

	public static CompoundGradient makeDataModelCompoundGradient(List<DataModel> dataModels) {
		List<GradientWrtParameterProvider> gpp = new ArrayList<>();
		for (DataModel dm : dataModels) {
			GradientProvider gp = (GradientProvider) dm.getLikelihood();
			gpp.add(new GradientWrtParameterProvider.ParameterWrapper(
					gp, dm.getCompoundParameter(), dm.getLikelihood()));
		}
		return new CompoundDerivative(gpp);
	}

	public static void main(String[] args) {

		int chainLength = 1100000;
		int burnIn = 100000;
		int subSampleFrequency = 1000;


		List<DataModel> allDataModels = new ArrayList<>();
		allDataModels.add(new ExtendingEmpiricalDataModel("ForDavid/grids_example_1.csv"));
		allDataModels.add(new ExtendingEmpiricalDataModel("ForDavid/grids_example_2.csv"));
		allDataModels.add(new ExtendingEmpiricalDataModel("ForDavid/grids_example_3.csv"));

		HierarchicalMetaAnalysisConfiguration cg = new HierarchicalMetaAnalysisConfiguration();
		//cg.exposureHyperStdDev = 0.0001; // fix exposure effect
		//cg.exposureHyperLocation = 0.5;

		//cg.hierarchicalLocationSecondaryHyperStdDev = 0.0001; // fix source.mean  to default 0

		//cg.tauShape = 1.0;
		//cg.tauScale = 100.0; // change up prior for tau, precision for the iid normal error term
		//cg.startingTau = 0.5;

		//cg.gammaPrimaryHyperShape = 1000000; // change up prior for across-outcome / across-datasource precision term
		//cg.gammaPrimaryHyperScale = 0.0001;

		cg.separateEffectPrior = true; // try with separate prior on main effect

		HierarchicalMetaAnalysis analysis = new HierarchicalMetaAnalysis(allDataModels,
				cg);

		Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, cg.seed);

		runner.run();

		runner.processSamples();
	}
}
