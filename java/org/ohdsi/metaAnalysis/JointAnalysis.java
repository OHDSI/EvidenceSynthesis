/*******************************************************************************
 * Copyright 2021 Observational Health Data Sciences and Informatics
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
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.inference.operators.*;
import dr.math.distributions.NormalDistribution;
import org.ohdsi.data.CoxData;
import org.ohdsi.data.SortedCoxData;
import org.ohdsi.likelihood.CoxPartialLikelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class JointAnalysis implements Analysis {

	private final Likelihood likelihood;
	private final Likelihood prior;
	private final Likelihood joint;
	private final OperatorSchedule schedule;

	private final Parameter beta;

    @SuppressWarnings("WeakerAccess")
	public JointAnalysis(CoxData data, double betaPriorSd) { this (data.getSortedData(), betaPriorSd); }

	@SuppressWarnings("WeakerAccess")
	public JointAnalysis(SortedCoxData data, double betaPriorSd) {

		// Build likelihood
		beta = new Parameter.Default("beta", 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
		likelihood = new CoxPartialLikelihood(beta, data);

		// Build prior
		DistributionLikelihood betaPrior = new DistributionLikelihood(new NormalDistribution(0, betaPriorSd));
		betaPrior.addData(beta);
		prior = betaPrior;

		// Build joint
		joint = new CompoundLikelihood(Arrays.asList(likelihood, prior));

		// Build transition kernel
		schedule = new SimpleOperatorSchedule(1000, 0.0);

		double defaultWeight = 1.0;
		AdaptationMode mode = AdaptationMode.ADAPTATION_ON;
		RandomWalkOperator.BoundaryCondition condition = RandomWalkOperator.BoundaryCondition.reflecting;

		schedule.addOperator(new RandomWalkOperator(beta, null, 0.75, condition, defaultWeight, mode));
	}

	@Override
	 public List<Loggable> getLoggerColumns() {

	    likelihood.setId("likelihood");
	    prior.setId("prior");

		List<Loggable> columns = new ArrayList<>();
		columns.add(likelihood);
		columns.add(prior);
		columns.add(beta);

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

	@SuppressWarnings("unused")
	private static SortedCoxData parseDataFile(String fileName) {
		return null;
	}

	public static void main(String[] args) {

		int chainLength = 1100000;
		int burnIn = 100000;
		int subSampleFrequency = 1000;

        // With strata
        int[] outcome = new int[] { 1, 1, 0, 1, 1, 0, 1 };
        double[] time = new double[] { 4, 3, 3, 2, 2, 1, 1 };
        double[] x = new double[] { 0, 2, 0, 0, 1, 1, 1 };
        int[] strata = new int[] { 0, 0, 1, 1, 1, 0, 0 };

        CoxData data = new CoxData(strata, outcome, time, x);

		JointAnalysis analysis = new JointAnalysis(data, 1000);

		Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, 666);

		runner.run();

		runner.processSamples();
	}
}
