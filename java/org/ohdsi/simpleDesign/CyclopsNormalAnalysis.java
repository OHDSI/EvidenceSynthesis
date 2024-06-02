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
import dr.inference.regression.CyclopsRegressionModel;
import dr.math.distributions.NormalDistribution;
import org.ohdsi.mcmc.Analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class CyclopsNormalAnalysis implements Analysis {

    private final Likelihood likelihood;
    private final Likelihood prior;
    private final Likelihood joint;

    private final Parameter beta;

    private final OperatorSchedule schedule;

    @SuppressWarnings("unused")
    public CyclopsNormalAnalysis(CyclopsRegressionModel cyclops, double priorMean, double priorSd,
                                 double startingValue) {

        // Build likelihood
        int threadCount = 1;
        likelihood = new CompoundLikelihood(threadCount, Collections.singleton(cyclops));
        beta = cyclops.getParameter();
        for (int i = 0; i < beta.getDimension(); ++i) {
            beta.setParameterValue(i, startingValue);
        }

        // Build prior
        DistributionLikelihood betaPrior = new DistributionLikelihood(new NormalDistribution(priorMean, priorSd));
        betaPrior.addData(beta);
        prior = betaPrior;
        prior.setId("prior");

        // Build joint
        joint = new CompoundLikelihood(Arrays.asList(likelihood, prior));
        joint.setId("joint");

        // Build transition kernel
        schedule = new SimpleOperatorSchedule(1000, 0.0);
        double defaultWeight = 1.0;
        AdaptationMode mode = AdaptationMode.ADAPTATION_ON;

        RandomWalkOperator.BoundaryCondition condition = RandomWalkOperator.BoundaryCondition.reflecting;
        schedule.addOperator(new RandomWalkOperator(beta, null, 0.75, condition, defaultWeight, mode));
    }

    @Override
    public List<Loggable> getLoggerColumns() {

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

}



