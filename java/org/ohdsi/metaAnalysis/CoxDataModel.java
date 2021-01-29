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

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import org.ohdsi.data.CoxData;
import org.ohdsi.likelihood.CoxPartialLikelihood;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class CoxDataModel implements DataModel {

    private final List<Parameter> thetaList;
    private final List<Likelihood> likelihoodList;

    private CompoundLikelihood likelihood = null;
    private Parameter theta = null;

    public CoxDataModel() {
        thetaList = new ArrayList<>();
        likelihoodList = new ArrayList<>();
   	}

    public void addLikelihoodData(int[] id,
                                  int[] outcome,
                                  double[] time,
                                  double[] covariate) {

        Parameter theta = new Parameter.Default("theta" + (thetaList.size() + 1), 0.1,
                        Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        CoxData data = new CoxData(id, outcome, time, covariate);
        CoxPartialLikelihood thisLikelihood = new CoxPartialLikelihood(theta, data.getSortedData());

        thetaList.add(theta);
        likelihoodList.add(thisLikelihood);
    }

    public void finish() {
        likelihood = new CompoundLikelihood(likelihoodList);
        theta = new CompoundParameter("theta", thetaList.toArray(new Parameter[] {}));
    }

    @Override
    public Likelihood getLikelihood() { return likelihood; }

    @Override
    public Parameter getCompoundParameter() { return theta; }

    @Override
    public List<Parameter> getIndividualParameters() { return thetaList; }
}
