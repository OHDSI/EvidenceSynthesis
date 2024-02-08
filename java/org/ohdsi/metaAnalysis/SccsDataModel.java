/*******************************************************************************
 * Copyright 2024 Observational Health Data Sciences and Informatics
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

import org.ohdsi.data.SccsData;
import org.ohdsi.likelihood.SccsPartialLikelihood;
import org.ohdsi.mcmc.Runner;

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

/**
 * @author Martijn Schuemie
 */
public class SccsDataModel implements DataModel {

    private final List<Parameter> thetaList;
    private final List<Likelihood> likelihoodList;
    private final List<Integer> identifierList;

    private CompoundLikelihood likelihood = null;
    private Parameter theta = null;

    public SccsDataModel() {
        thetaList = new ArrayList<>();
        likelihoodList = new ArrayList<>();
        identifierList = new ArrayList<>();
   	}

    public void addLikelihoodData(int[] y, double[] a, double[][] x, int[] stratumId, double[] time) {

        Parameter theta = new Parameter.Default("theta" + (thetaList.size() + 1), 0.1,
                        Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        SccsData data = new SccsData(y, a, x, stratumId, time);
        SccsPartialLikelihood thisLikelihood = new SccsPartialLikelihood(theta, data);

        thetaList.add(theta);
        likelihoodList.add(thisLikelihood);
        identifierList.add(identifierList.size() + 1);
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

    @Override
    public List<Integer> getIdentifiers() { return identifierList; }
    
    public static void main(String[] args) {
    	SccsDataModel dataModel = new SccsDataModel();
    	int[] y = new int[] { 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0 };
		double[] time = new double[] { 107, 21, 54, 183, 41, 21, 120, 183, 78, 21, 83, 183, 82, 21, 79, 183, 81, 21, 80, 183, 44, 21, 117, 183, 119, 21, 42, 183, 145, 21, 16, 183, 77, 21, 84, 183, 182, 183 };
		double[] a = new double[] { 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0 };
		double[][] x = new double[][] { {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {1}, {0}, {1}};
		int[] stratumId = new int[] { 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10 };
    	dataModel.addLikelihoodData(y, a, x, stratumId, time);
//    	dataModel.addLikelihoodData(y, a, x, stratumId, time);
    	dataModel.finish();
    	HalfNormalOnStdDevPrior prior = new HalfNormalOnStdDevPrior(0, 0.5);
    	MetaAnalysis metaAnalysis = new MetaAnalysis(dataModel, prior, 2.0);
    	Runner runner = new Runner(metaAnalysis, 10, 1, 1, 123, false);
    	runner.run();
    }
}
