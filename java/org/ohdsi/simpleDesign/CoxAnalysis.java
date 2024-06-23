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
package org.ohdsi.simpleDesign;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.RandomWalkOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.math.distributions.NormalDistribution;
import org.ohdsi.data.CoxData;
import org.ohdsi.data.SortedCoxData;
import org.ohdsi.likelihood.MultivariableCoxPartialLikelihood;
import org.ohdsi.mcmc.Analysis;
import org.ohdsi.mcmc.Runner;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.ohdsi.likelihood.MultivariableCoxPartialLikelihood.exampleMultivariableData;

public class CoxAnalysis implements Analysis {

    private final Likelihood likelihood;
    private final Likelihood prior;
    private final Likelihood joint;
    private final OperatorSchedule schedule;

    private final Parameter beta;

    public CoxAnalysis(CoxData data, double betaPriorMean, double betaPriorSd) {
        this(data.getSortedData(), betaPriorMean, betaPriorSd);
    }

    public CoxAnalysis(SortedCoxData data, double betaPriorMean, double betaPriorSd) {

        // Build likelihood
        beta = Analysis.makeParameter("beta", data.getCovariateDimension());
        likelihood = new MultivariableCoxPartialLikelihood(beta, data);

        // Build prior
        DistributionLikelihood betaPrior = new DistributionLikelihood(new NormalDistribution(betaPriorMean, betaPriorSd));
        betaPrior.addData(beta);

        prior = new CompoundLikelihood(List.of(betaPrior));

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

    public static void main(String[] args) {

        int chainLength = 11000;
        int burnIn = 1000;
        int subSampleFrequency = 10;

        SortedCoxData data = exampleMultivariableData();

        CoxAnalysis analysis = new CoxAnalysis(data, 0, 1);

        Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, 666);

        runner.run();

        System.out.println();
        runner.processSamples();
    }
}
