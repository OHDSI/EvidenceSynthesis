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

import dr.inference.model.*;
import dr.math.distributions.NormalDistribution;

import java.util.List;

public class SimpleLinearModel extends AbstractModelLikelihood {

    private final Parameter argument;
    private final DesignMatrix designMatrix;
    private final Parameter effects;
    private final Parameter precision;
    private final List<Integer> betaToTauIndexMap; // For heteroscedastic model

    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;

    private boolean innerProductKnown;
    private boolean storedInnerProductKnown;

    private double logLikelihood;
    private double storedLogLikelihood;

    private double[] innerProduct;
    private double[] storedInnerProduct;

    /**
     * Original constructor for a homoscedastic linear model.
     */
    public SimpleLinearModel(String name,
                             Parameter argument,
                             DesignMatrix designMatrix,
                             Parameter effects,
                             Parameter precision) {
        // Call the main constructor with a null map for backward compatibility
        this(name, argument, designMatrix, effects, precision, null);
    }

    /**
     * Main constructor supporting both homoscedastic and heteroscedastic models.
     *
     * @param name               The model name.
     * @param argument           The parameter representing the data points (e.g., betas).
     * @param designMatrix       The design matrix.
     * @param effects            The parameter representing the effects.
     * @param precision          A parameter (or compound parameter) of precisions.
     * @param betaToTauIndexMap  A map from each argument to a precision index. If null, model is homoscedastic.
     */
    public SimpleLinearModel(String name,
                             Parameter argument,
                             DesignMatrix designMatrix,
                             Parameter effects,
                             Parameter precision,
                             List<Integer> betaToTauIndexMap) {
        super(name);

        this.argument = argument;
        this.designMatrix = designMatrix;
        this.effects = effects;
        this.precision = precision;
        this.betaToTauIndexMap = betaToTauIndexMap;

        if (designMatrix.getRowDimension() != argument.getDimension()) {
            throw new IllegalArgumentException("Invalid parameter dimensions: design matrix rows must match argument dimension");
        }

        if (designMatrix.getColumnDimension() != effects.getDimension()) {
            throw new IllegalArgumentException("Invalid parameter dimensions: design matrix columns must match effects dimension");
        }

        // If the map is provided, validate its size
        if (betaToTauIndexMap != null && betaToTauIndexMap.size() != argument.getDimension()) {
            throw new IllegalArgumentException("Invalid parameter dimensions: betaToTauIndexMap size must match argument dimension");
        }

        addVariable(argument);
        addVariable(designMatrix);
        addVariable(effects);
        addVariable(precision);

        innerProduct = new double[designMatrix.getRowDimension()];
        storedInnerProduct = new double[designMatrix.getRowDimension()];

        likelihoodKnown = false;
        innerProductKnown = false;
    }

    private void computeInnerProduct(double[] product) {
        for (int i = 0; i < designMatrix.getRowDimension(); ++i) {
            product[i] = 0.0;
            for (int j = 0; j < designMatrix.getColumnDimension(); ++j) {
                product[i] += designMatrix.getParameterValue(i, j) * effects.getParameterValue(j);
            }
        }
    }

    private double calculateLogLikelihood() {

        if (!innerProductKnown) {
            computeInnerProduct(innerProduct);
            innerProductKnown = true;
        }

        double logLikelihood = 0.0;

        if (betaToTauIndexMap != null) {
            // Heteroscedastic case: Use the map to find the correct precision for each argument
            for (int i = 0; i < argument.getDimension(); ++i) {
                int precisionIndex = betaToTauIndexMap.get(i);
                double tau = precision.getParameterValue(precisionIndex);
                double sd = 1.0 / Math.sqrt(tau);
                logLikelihood += NormalDistribution.logPdf(argument.getParameterValue(i), innerProduct[i], sd);
            }
        } else {
            // Homoscedastic case: Use a single precision for all arguments
            double tau = precision.getParameterValue(0);
            double sd = 1.0 / Math.sqrt(tau);
            for (int i = 0; i < argument.getDimension(); ++i) {
                logLikelihood += NormalDistribution.logPdf(argument.getParameterValue(i), innerProduct[i], sd);
            }
        }
        return logLikelihood;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object o, int i) {
        throw new RuntimeException("Should not occur");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int i, Variable.ChangeType changeType) {
        if (variable == precision || variable == argument) {
            likelihoodKnown = false;
        } else if (variable == designMatrix || variable == effects) {
            likelihoodKnown = false;
            innerProductKnown = false;
        } else {
            throw new RuntimeException("Should not get here");
        }
    }

    @Override
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;

        storedInnerProductKnown = innerProductKnown;

        if (innerProductKnown) {
            System.arraycopy(innerProduct, 0, storedInnerProduct, 0, innerProduct.length);
        }
    }

    @Override
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;

        innerProductKnown = storedInnerProductKnown;

        if (innerProductKnown) {
            double[] swap = innerProduct;
            innerProduct = storedInnerProduct;
            storedInnerProduct = swap;
        }
    }

    @Override
    protected void acceptState() {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        innerProductKnown = false;
    }

    public final Parameter getArgument() { return argument; }

    public final DesignMatrix getDesignMatrix() { return designMatrix; }

    public final Parameter getEffects() { return effects; }

    public final Parameter getPrecision() { return precision; }

    public final double[] getInnerProduct() {
        if (!innerProductKnown) {
            computeInnerProduct(innerProduct);
            innerProductKnown = true;
        }
        return innerProduct;
    }

    public final List<Integer> getBetaToTauIndexMap() {
        return betaToTauIndexMap;
    }
}
