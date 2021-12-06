package org.ohdsi.simpleDesign;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.inference.operators.*;
import dr.math.distributions.NormalDistribution;
import org.ohdsi.likelihood.EmpiricalModelLikelihood;
import org.ohdsi.mcmc.Analysis;
import org.ohdsi.mcmc.Runner;
import org.ohdsi.metaAnalysis.EmpiricalDataModel;
import org.ohdsi.metaAnalysis.ExtendingEmpiricalDataModel;

@SuppressWarnings("unused")
public class ProfileNormalAnalysis implements Analysis {

    private final Likelihood likelihood;
    private final Likelihood prior;
    private final Likelihood joint;

    private final Parameter beta;

    private final OperatorSchedule schedule;

    @SuppressWarnings("unused")
    public ProfileNormalAnalysis(EmpiricalDataModel profileLikelihood, double priorMean, double priorSd,
                                 double startingValue) {

        // Build likelihood
        likelihood = new EmpiricalModelLikelihood("likelihood", profileLikelihood);
        beta = profileLikelihood.getCompoundParameter();
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
        for (Parameter p : profileLikelihood.getIndividualParameters()) {
            schedule.addOperator(new RandomWalkOperator(p, null, 0.75, condition, defaultWeight, mode));
        }
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

    public static void main(String[] args) {

        int chainLength = 1100000;
        int burnIn = 100000;
        int subSampleFrequency = 1000;

        EmpiricalDataModel dataModel = new ExtendingEmpiricalDataModel("profile.txt");

        Analysis analysis = new ProfileNormalAnalysis(dataModel, 0, 10, 0);

        Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, 666);
        runner.run();
        runner.processSamples();
    }
}



