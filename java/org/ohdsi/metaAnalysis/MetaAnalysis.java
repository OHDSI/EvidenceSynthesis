/*******************************************************************************
 * Copyright 2020 Observational Health Data Sciences and Informatics
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.distributions.NormalDistribution;

public class MetaAnalysis implements Analysis {

	private final Likelihood likelihood;
	private final Likelihood prior;
	private final Likelihood joint;

	private final Parameter theta;
	private final Parameter mu;
	private final Parameter tau;

	private final OperatorSchedule schedule;

	public MetaAnalysis(DataModel dataModel, ScalePrior scalePrior, double muPriorSd) {

		// Build likelihood
		Likelihood dataLikelihood = dataModel.getLikelihood();
		theta = dataModel.getCompoundParameter();

		mu = new Parameter.Default("mu", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
		tau = scalePrior.getParameter();
		boolean isPrecision = scalePrior.isPrecision();

		DistributionLikelihood hierarchicalLikelihood = new DistributionLikelihood(
				new NormalDistributionModel(mu, tau, isPrecision));
		hierarchicalLikelihood.addData(theta);

		int defaultThreads = 0; // No thread pools
		likelihood = new CompoundLikelihood(defaultThreads, Arrays.asList(dataLikelihood, hierarchicalLikelihood));

		// Build prior
		DistributionLikelihood muPrior = new DistributionLikelihood(new NormalDistribution(0, muPriorSd));
		muPrior.addData(mu);

		Likelihood tauPrior = scalePrior.getPrior();

		prior = new CompoundLikelihood(Arrays.asList(muPrior, tauPrior));

		// Build joint
		joint = new CompoundLikelihood(Arrays.asList(likelihood, prior));

		// Build transition kernel
		schedule = new SimpleOperatorSchedule(1000, 0.0);
		double defaultWeight = 1.0;
		schedule.addOperator(
				new NormalNormalMeanGibbsOperator(hierarchicalLikelihood, muPrior.getDistribution(), defaultWeight));

		AdaptationMode mode = AdaptationMode.ADAPTATION_ON;

		schedule.addOperator(scalePrior.getOperator(hierarchicalLikelihood, defaultWeight, mode));

		RandomWalkOperator.BoundaryCondition condition = RandomWalkOperator.BoundaryCondition.reflecting;
		for (Parameter p : dataModel.getIndividualParameters()) {
			schedule.addOperator(new RandomWalkOperator(p, null, 0.75, condition, defaultWeight, mode));
		}
	}

	@Override
	 public List<Loggable> getLoggerColumns() {

		List<Loggable> columns = new ArrayList<>();
		columns.add(likelihood);
		columns.add(prior);
		columns.add(mu);
		columns.add(tau);
		if (theta != null) {
			columns.add(theta);
		}

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

	public static void main(String[] args) {

		int chainLength = 1100000;
		int burnIn = 100000;
		int subSampleFrequency = 1000;

//        MetaAnalysis analysis = new MetaAnalysis(
//                new DataModel.Empirical("c:/temp/simGridData.csv"),
//                new ScalePrior.GammaOnPrecision(0.001, 1000.0),
//                chainLength, burnIn, subSampleFrequency);

//		MetaAnalysis analysis = new MetaAnalysis(new ExtendingEmpiricalDataModel("c:/temp/grids_example_3.csv"),
//				new HalfNormalOnStdDevPrior(0.0, 2), 1000, chainLength, burnIn, subSampleFrequency);

//		MetaAnalysis analysis = new MetaAnalysis(new NormalDataModel("c:/temp/normal_example_1.csv"),
//				new HalfCauchyOnStdDevPrior(0.0, 2), chainLength, burnIn, subSampleFrequency);

		MetaAnalysis analysis = new MetaAnalysis(new SkewNormalDataModel("c:/temp/skewnormal_example_3.csv"),
				new HalfNormalOnStdDevPrior(0.0, 2), 1000);

		Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency);

		runner.run();

		runner.processSamples();
	}
}
