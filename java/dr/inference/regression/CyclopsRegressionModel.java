/*
 * CyclopsRegressionModel.java
 *
 * Copyright (c) 2002-2024 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.regression;

import dr.inference.model.*;
import org.ohdsi.metaAnalysis.DataModel;

import java.util.*;

import static org.ohdsi.mcmc.Analysis.makeParameter;

// TODO: Do we want to include the Cyclops regularization? (may be helpful for `findMode()`)

/**
 * @author Marc A. Suchard
 */
public class CyclopsRegressionModel extends AbstractModelLikelihood implements DataModel {

    private final RegressionInCyclops cyclops;
    private final Parameter beta;
    private final int cyclopsDim;
    private final int cyclopsOffset;

    private double logLikelihood;
    private double storedLogLikelihood;

    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;

    private boolean betaKnown;
    //    private boolean storedBetaKnown;
    private boolean updateAllBeta;

    private final Queue<Integer> betaDimChanged = new LinkedList<>();

    public CyclopsRegressionModel(String name, String libraryFileName, int instance,
                                  int dim, boolean useCyclopsStartingValue) {
        this(name, libraryFileName, instance, makeParameter(name + ".beta", dim), useCyclopsStartingValue);
    }

    public CyclopsRegressionModel(String name, String libraryFileName, int instance,
                                  Parameter beta, boolean useCyclopsStartingValues) {
        super(name);
        this.cyclops = new RegressionInCyclops(libraryFileName, instance);

        this.beta = beta;
        this.cyclopsOffset = cyclops.hasOffset();
        this.cyclopsDim = cyclops.getBetaSize();
        if (beta.getDimension() != cyclopsDim - cyclopsOffset) {
            throw new IllegalArgumentException("Invalid beta parameter");
        }

        if (useCyclopsStartingValues) {
            setBetaParameterFromCyclops();
            betaKnown = true;
        } else {
            betaKnown = false;
        }

        this.updateAllBeta = false;

        addVariable(beta);

        beta.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
                beta.getDimension()));
    }

    @SuppressWarnings("unused")
    public void findMode() {
        cyclops.findMode();
        setBetaParameterFromCyclops();

        betaKnown = true;
        likelihoodKnown = false;
    }

    private void setBetaParameterFromCyclops() {
        double[] values = new double[cyclopsDim];

        cyclops.getBeta(values);

        for (int i = 0; i < cyclopsDim - cyclopsOffset; ++i) {
            beta.setParameterValueQuietly(i, values[i + cyclopsOffset]);
        }
        beta.fireParameterChangedEvent();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object o, int i) {
        throw new IllegalArgumentException("Unknown model");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType changeType) {
        if (variable == beta) {
            betaKnown = false;
            likelihoodKnown = false;
            if (changeType == Variable.ChangeType.ALL_VALUES_CHANGED) {
                updateAllBeta = true;
                betaDimChanged.clear();
            } else {
                betaDimChanged.add(index);
            }
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    @Override
    protected void storeState() {
//        storedBetaKnown = betaKnown;
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
//        betaKnown = storedBetaKnown;
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    @Override
    protected void acceptState() {
    }

    @Override
    public Model getModel() {
        return this;
    }

    public Parameter getParameter() {
        return beta;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    private double calculateLogLikelihood() {
        if (!betaKnown) {
            setBetaInCyclops();
            betaKnown = true;
            updateAllBeta = false;
        }

        return cyclops.getLogLikelihood();
    }

    private void setBetaInCyclops() {
        if (updateAllBeta || betaDimChanged.isEmpty()) {
            double[] values = beta.getParameterValues();

            if (cyclopsOffset > 0) {
                double[] augmented = new double[cyclopsDim];
                for (int i = 0; i < cyclopsOffset; ++i) {
                    augmented[i] = cyclops.getBeta(i);
                }
                System.arraycopy(values, 0, augmented, cyclopsOffset, cyclopsDim - cyclopsOffset);
                values = augmented;
            }

            cyclops.setBeta(values);
            betaDimChanged.clear();
        } else {
            while (!betaDimChanged.isEmpty()) {
                final int index = betaDimChanged.remove();
                cyclops.setBeta(index + cyclopsOffset, beta.getParameterValue(index));
            }
        }
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        betaKnown = false;
        cyclops.makeDirty();
    }

    protected double[] getGradientWrtBeta() {
        double[] gradient = new double[cyclopsDim];
        cyclops.getLogLikelihoodGradient(gradient);
        return gradient;
    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Parameter getCompoundParameter() {
        return beta;
    }

    @Override
    public List<Parameter> getIndividualParameters() {
        return Collections.singletonList(beta);
    }

    @Override
    public List<Integer> getIdentifiers() {
        return null;
    }
}
