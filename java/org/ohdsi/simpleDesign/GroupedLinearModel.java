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

import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.model.*;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;

public class GroupedLinearModel extends SimpleLinearModel implements ConjugateWishartStatisticsProvider {

    public enum Grouping {
        BY_ROW {
            @Override
            public int getDimension(MatrixParameterInterface argument) {
                return argument.getColumnDimension();
            }

            @Override
            public int index(int group, int dimensionWithinGroup, int numGroups, int mvnDim) {
                return dimensionWithinGroup * numGroups + group;
            }
        },
        BY_COLUMN {
            @Override
            public int getDimension(MatrixParameterInterface argument) {
                return argument.getRowDimension();
            }

            @Override
            public int index(int group, int dimensionWithinGroup, int numGroups, int mvnDim) {
                return group * mvnDim + dimensionWithinGroup;
            }
        };

        public abstract int getDimension(MatrixParameterInterface argument);

        public abstract int index(int group, int dimensionWithinGroup, int numGroups, int mvnDim);
    }

    private final Grouping grouping;
    private final int numGroups;
    private final int mvnDim;
    private final MatrixParameterInterface precision;
    private final MultivariateNormalDistributionModel mvn;

    private final double[] delta;

    public GroupedLinearModel(String name,
                              MatrixParameter argument,
                              DesignMatrix designMatrix,
                              Parameter effects,
                              MatrixParameter precision,
                              Grouping grouping) {

        super(name, argument, designMatrix, effects, precision);
        this.mvnDim = grouping.getDimension(argument);

        assert mvnDim == precision.getRowDimension();
        assert mvnDim == precision.getColumnDimension();

        this.mvn = new MultivariateNormalDistributionModel(
                new Parameter.Default(mvnDim, 0.0), precision);

        this.precision = precision;
        this.grouping = grouping;
        this.numGroups = argument.getDimension() / mvnDim;
        this.delta = new double[mvnDim];

        addModel(mvn);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object o, int i) {
        if (model == mvn) {
            likelihoodKnown = false;
        } else {
            throw new RuntimeException("Should not occur");
        }
    }

    @Override
    protected double calculateLogLikelihood() {

        checkInnerProduct();

        double logLikelihood = 0.0;
        for (int g = 0; g < numGroups; ++g) {
            logLikelihood += mvn.logPdf(getDeltaForGroup(g));
        }

        return logLikelihood;
    }

    private double[] getDeltaForGroup(int group) {

        for (int j = 0; j < mvnDim; ++j) {
            int index = grouping.index(group, j, numGroups, mvnDim);
            delta[j] = argument.getParameterValue(index) - innerProduct[index];
        }

        return delta;
    }

    private int getDf() { return numGroups; }

    private double[] getOuterProducts() {

        checkInnerProduct();

        double[] outerProducts = new double[mvnDim * mvnDim];

        for (int g = 0; g < numGroups; ++g) {

            double[] delta = getDeltaForGroup(g);
            int index = 0;
            for (int i = 0; i < mvnDim; ++i) {
                for (int j = 0; j < mvnDim; ++j) {
                    outerProducts[index] += delta[i] * delta[j];
                    ++index;
                }
            }
        }

        return outerProducts;
    }

    @Override
    public WishartSufficientStatistics getWishartStatistics() {
        int df = getDf();
        double[] outerProducts = getOuterProducts();
        return new WishartSufficientStatistics(df, outerProducts);
    }

    @Override
    public MatrixParameterInterface getPrecisionParameter() {
        return precision;
    }
}
