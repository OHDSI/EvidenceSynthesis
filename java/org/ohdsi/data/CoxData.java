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
package org.ohdsi.data;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import org.ohdsi.likelihood.CoxPartialLikelihood;

import java.util.*;

/**
 * @author Marc A. Suchard
 */
public class CoxData {

    private final SortedCoxData data;

    public CoxData(int[] outcome, double[] time, double[] covariate) {
        this(null, outcome, time, covariate);
    }

    public CoxData(int[] id, int[] outcome, double[] time, double[] covariate) {

        if (outcome.length != time.length || outcome.length != covariate.length ||
                (id != null && id.length != outcome.length)) {
            throw new IllegalArgumentException("All dimensions must be equal");
        }

        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < outcome.length; ++i) {
            indices.add(i);
        }

        indices.sort(Comparator.comparing(o -> (outcome[o])));
        indices.sort(Comparator.comparing(o -> (-time[o])));
        if (id != null) {
            indices.sort(Comparator.comparing(o -> (id[o])));
        }

        int[] y = getInt(indices, outcome);
        double[] x = getDouble(indices, covariate);
        int[] p = getInt(indices, id);
        double[] t = getDouble(indices, time);
        double[] w = getWeights(y, p, t);

        this.data = new SortedCoxData(y, x, getStrata(p, y.length), w);
    }

    public SortedCoxData getSortedData() { return data; }

    private boolean isFailureTie(int[] y, double[] t, int[] id, int i) {
        boolean failureTie = (y[i] == 1 && y[i + 1] == 1 && t[i] == t[i + 1]);
        if (id != null) {
            failureTie = failureTie && id[i] == id[i + 1];
        }
        return failureTie;
    }

    private double[] getWeights(int[] y, int[] id, double[] t) {
        double[] weights = new double[y.length];

        double w = 1.0;
        for (int i = 0; i < y.length - 1; ++i) {
            if (isFailureTie(y, t, id, i)) {
                w += 1.0; // Breslow approximation for failure ties
                weights[i] = 0.0;
            } else {
                weights[i] = y[i] == 1 ? w : 0.0;
                w = 1.0;
            }
        }
        weights[y.length - 1] = y[y.length - 1] == 1 ? w : 0.0;

        return weights;
    }

    private int[] getInt(List<Integer> indices, int[] outcome) {
        if (outcome == null) {
            return null;
        }
        
        int[] y = new int[outcome.length];
        for (int i = 0; i < y.length; ++i) {
            y[i] = outcome[indices.get(i)];
        }
        return y;
    }

    private double[] getDouble(List<Integer> indices, double[] covariate) {
        double[] x = new double[covariate.length];
        for (int i = 0; i < x.length; ++i) {
            x[i] = covariate[indices.get(i)];
        }
        return x;
    }

    private int[] getStrata(int[] id, int length) {
        if (id == null) {
            return new int[] { length };
        }

        List<Integer> strata = new ArrayList<>();
        int last = id[0];
        int i = 1;
        for ( ; i < id.length; ++i) {
            int current = id[i];
            if (current != last) {
                strata.add(i);
                last = current;
            }
        }
        strata.add(i);

        return strata.stream().mapToInt(index -> index).toArray();
    }
    
    public static void main(String[] args) {

        // No strata
//   		int[] outcome = new int[] { 1, 1, 0, 1, 1, 0, 1 };
//   		double[] time = new double[] { 4, 3.5, 3, 2.5, 2, 1.5, 1 };
//   		double[] covariate = new double[] { 0, 2, 0, 0, 1, 1, 1 };
////   		int[] strata = new int[] { 7 };
//   		double beta = 0.3883064;

        // With strata
//        int[] y = new int[] { 1, 1, 0, 1, 0, 1, 1 };
//        double[] x = new double[] { 0, 2, 1, 1, 0, 0, 1 };
//        int[] strata = new int[] { 4, 7 };
//        double beta = 1.205852;

        int[] outcome = new int[] { 1, 1, 0, 1, 1, 0, 1 };
        double[] time = new double[] { 4, 3, 3, 2, 2, 1, 1 };
        double[] covariate = new double[] { 0, 2, 0, 0, 1, 1, 1 };
        int[] id = new int[] { 0, 0, 1, 1, 1, 0, 0 };
        double beta = 0.7357498;
        // log-Lik = -3.805248

   		CoxData data = new CoxData(id, outcome, time, covariate);

   		Parameter parameter = new Parameter.Default(beta);

   		Likelihood cox = new CoxPartialLikelihood(parameter, data.getSortedData());

   		System.err.println(cox.getLogLikelihood());
   	}
}
