package org.ohdsi.metaAnalysis;

import java.util.ArrayList;
import java.util.List;

import dr.inference.distribution.EmpiricalDistributionData;
import dr.inference.distribution.EmpiricalDistributionLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

public interface DataModel {

    Likelihood getLikelihood();

    Parameter getCompoundParameter();

    List<Parameter> getIndividualParameters();

    abstract class Base implements DataModel {

        private final List<Parameter> thetaList = new ArrayList<Parameter>();
        private EmpiricalDistributionLikelihood likelihood = null;
        private Parameter theta = null;
        private List<EmpiricalDistributionData> dataList = new ArrayList<EmpiricalDistributionData>();

        public void addLikelihoodParameters(double[] x, double[] ll) {
            dataList.add(new EmpiricalDistributionData(x, ll, true));
            thetaList.add(new Parameter.Default("theta" + (thetaList.size() + 1), 0.1,
                    Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY));
        }

        public void finish() {
            likelihood = makeFunctionalForm(dataList);
            theta = new CompoundParameter("theta", thetaList.toArray(new Parameter[] {}));
            likelihood.addData(theta);
        }

        abstract EmpiricalDistributionLikelihood makeFunctionalForm(List<EmpiricalDistributionData> dataList);

        @Override
        public Likelihood getLikelihood() { return likelihood; }

        @Override
        public Parameter getCompoundParameter() { return theta; }

        @Override
        public List<Parameter> getIndividualParameters() { return thetaList; }
    }
}
