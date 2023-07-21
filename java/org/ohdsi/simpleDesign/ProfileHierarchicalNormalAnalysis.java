package org.ohdsi.simpleDesign;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.RandomWalkOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.NormalDistribution;
import org.ohdsi.likelihood.CachedModelLikelihood;
import org.ohdsi.mcmc.Analysis;
import org.ohdsi.mcmc.Runner;
import org.ohdsi.metaAnalysis.EmpiricalDataModel;
import org.ohdsi.metaAnalysis.ExtendingEmpiricalDataModel;

import java.util.*;

@SuppressWarnings("unused")
public class ProfileHierarchicalNormalAnalysis implements Analysis {

    private final Likelihood likelihood;
    private final Likelihood prior;
    private final Likelihood joint;

    private final Parameter beta;
    private final Parameter mean;
    private final Parameter precision;

    private final OperatorSchedule schedule;


//    public ProfileHierarchicalNormalAnalysis(List<EmpiricalDataModel> profileLikelihoods, double priorEffectMean, double priorEffectSd,
//                                             double startingValue) {
//
//        // Build likelihood
//        this.likelihood = new EmpiricalModelLikelihood("likelihood", profileLikelihoods);
//        List<Likelihood> empiricalModels = new ArrayList<>();
//
//
//        this.beta = new CompoundParameter("beta");
//
//        for (EmpiricalDataModel dataModel : profileLikelihoods) {
//            empiricalModels.add(new EmpiricalModelLikelihood("likelihood", dataModel));
//            Parameter individualBeta = dataModel.getCompoundParameter();
//
//            if (individualBeta.getDimension() != 1) {
//                throw new RuntimeException("Not yet implemented");
//            }
//            individualBeta.setParameterValue(0, startingValue);
//
//            ((CompoundParameter) beta).addParameter(individualBeta);
//        }
//
//    }
    @SuppressWarnings("unused")
    public ProfileHierarchicalNormalAnalysis(List<EmpiricalDataModel> profileLikelihoods, double priorEffectMean, double priorEffectSd,
                                             double startingValue) {

        // Build likelihood
        List<Likelihood> empiricalModels = new ArrayList<>();
        this.beta = new CompoundParameter("beta");

        for (EmpiricalDataModel dataModel : profileLikelihoods) {
            empiricalModels.add(new CachedModelLikelihood("likelihood", dataModel));
            Parameter individualBeta = dataModel.getCompoundParameter();

            if (individualBeta.getDimension() != 1) {
                throw new RuntimeException("Not yet implemented");
            }
            individualBeta.setParameterValue(0, startingValue);

            ((CompoundParameter)beta).addParameter(individualBeta);
        }

        this.likelihood = new CompoundLikelihood(empiricalModels);

        this.mean = new Parameter.Default("mean", 1);
        mean.setParameterValue(0, 0.0);

        this.precision = new Parameter.Default("precision", 1);
        precision.setParameterValue(0, 1.0);

        DistributionLikelihood effectDistribution = new DistributionLikelihood(
                new NormalDistributionModel(mean, precision, true));
        effectDistribution.addData(beta);


        // Build prior
        DistributionLikelihood meanPrior = new DistributionLikelihood(new NormalDistribution(priorEffectMean, priorEffectSd));
        meanPrior.addData(mean);

        DistributionLikelihood precisionPrior = new DistributionLikelihood(new GammaDistribution(priorEffectMean, priorEffectSd));
        precisionPrior.addData(precision);

        this.prior = new CompoundLikelihood(Arrays.asList(effectDistribution, meanPrior, precisionPrior));
        prior.setId("prior");

        // Build joint
        joint = new CompoundLikelihood(Arrays.asList(likelihood, prior));
        joint.setId("joint");

        // Build transition kernel
        schedule = new SimpleOperatorSchedule(1000, 0.0);
        double defaultWeight = 1.0;
        AdaptationMode mode = AdaptationMode.ADAPTATION_ON;

        RandomWalkOperator.BoundaryCondition condition = RandomWalkOperator.BoundaryCondition.reflecting;

        for (EmpiricalDataModel dataModel : profileLikelihoods) {
            for (Parameter p : dataModel.getIndividualParameters()) {
                schedule.addOperator(new RandomWalkOperator(p, null, 0.75, condition, defaultWeight, mode));
            }
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

        Analysis analysis = new ProfileHierarchicalNormalAnalysis(Collections.singletonList(dataModel),
                0, 10, 0);

        Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, 666);
        runner.run();
        runner.processSamples();
    }
}



