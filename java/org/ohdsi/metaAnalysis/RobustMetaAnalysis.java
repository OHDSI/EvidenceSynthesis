/*******************************************************************************
 * Copyright 2022 Observational Health Data Sciences and Informatics
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
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.distribution.TDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.distributions.Distribution;
import org.ohdsi.mcmc.Runner;

public class RobustMetaAnalysis extends MetaAnalysis {


	public RobustMetaAnalysis(DataModel dataModel, ScalePrior scalePrior, double muPriorSd) {
		super(dataModel, scalePrior, muPriorSd);
	}

	protected MCMCOperator getMuOperator(Parameter mu,
										 DistributionLikelihood likelihood,
										 Distribution prior,
										 double weight) {
		AdaptationMode mode = AdaptationMode.ADAPTATION_ON;
		RandomWalkOperator.BoundaryCondition condition = RandomWalkOperator.BoundaryCondition.reflecting;

		return new RandomWalkOperator(mu, null, 0.75, condition, weight, mode);
	}

	protected ParametricDistributionModel getMuDistribution(Parameter mu, Parameter tau, boolean isPrecision) {
		return new TDistributionModel(mu, tau,
				new Parameter.Default("df", 4.0, 0.0, Double.POSITIVE_INFINITY));
	}

	public static void main(String[] args) {

		int chainLength = 1100000;
		int burnIn = 100000;
		int subSampleFrequency = 1000;

		RobustMetaAnalysis analysis = new RobustMetaAnalysis(new SkewNormalDataModel("c:/temp/skewnormal_example_3.csv"),
				new HalfNormalOnStdDevPrior(0.0, 2), 1000);

		Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, 666);

		runner.run();

		runner.processSamples();
	}
}
