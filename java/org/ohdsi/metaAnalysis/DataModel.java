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

    List<Integer> getIdentifiers();

    abstract class Base implements DataModel {

        private final List<Parameter> thetaList = new ArrayList<>();
        private EmpiricalDistributionLikelihood likelihood = null;
        private Parameter theta = null;
        private final List<EmpiricalDistributionData> dataList = new ArrayList<>();
        private final List<Integer> identifierList = new ArrayList<>();

        public void addLikelihoodParameters(double[] x, double[] ll) {
            addLikelihoodParameters(x, ll, identifierList.size() + 1);
        }

        public void addLikelihoodParameters(double[] x, double[] ll, int identifier) {
            dataList.add(new EmpiricalDistributionData(x, ll, true));
            thetaList.add(new Parameter.Default("theta" + (thetaList.size() + 1), 0.1,
                    Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY));
            identifierList.add(identifier);
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

        @Override
        public List<Integer> getIdentifiers() {  return identifierList; }
    }
}
