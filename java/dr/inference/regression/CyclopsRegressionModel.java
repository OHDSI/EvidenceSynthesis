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

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.WrappedVector;
import org.rosuda.JRI.*;

import java.util.LinkedList;
import java.util.Queue;

// TODO: Do we want to include the Cyclops regularization? (may be helpful for `findMode()`)

/**
 * @author Marc A. Suchard
 */
public class CyclopsRegressionModel extends AbstractModelLikelihood {

    private final RegressionInCyclops cyclops;
    private final Parameter beta;
    private final int dim;

    private double logLikelihood;
    private double storedLogLikelihood;

    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;

    private boolean betaKnown;
//    private boolean storedBetaKnown;
    private boolean updateAllBeta;

    private final Queue<Integer> betaDimChanged = new LinkedList<>();

    public CyclopsRegressionModel(String name, String libraryFileName, int instance,
                                  Parameter beta, boolean useCyclopsStartingValues) {
        super(name);
        this.cyclops = new RegressionInCyclops(libraryFileName, instance);

        this.beta = beta;
        this.dim = cyclops.getBetaSize();
        if (beta.getDimension() != dim) {
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
    }

    @SuppressWarnings("unused")
    public void findMode() {
        cyclops.findMode();
        setBetaParameterFromCyclops();

        betaKnown = true;
        likelihoodKnown = false;
    }

    private void setBetaParameterFromCyclops() {
        double[] values = new double[dim];
        cyclops.getBeta(values);

        for (int i = 0; i < dim; ++i) {
            beta.setParameterValueQuietly(i, values[i]);
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
    protected void acceptState() { }

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
            cyclops.setBeta(beta.getParameterValues());
        } else {
            while (!betaDimChanged.isEmpty()) {
                final int index = betaDimChanged.remove();
                cyclops.setBeta(index, beta.getParameterValue(index));
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
        double[] gradient = new double[dim];
        cyclops.getLogLikelihoodGradient(gradient);
        return gradient;
    }

    public static void main(String[] args) {

        Rengine rEngine = new Rengine(new String[]{"--no-save"}, false, null);

        rEngine.eval("library(Cyclops)");
        rEngine.eval("" +
                "dobson <- data.frame( " +
                "  counts = c(18,17,15,20,10,20,25,13,12)," +
                "  outcome = gl(3,1,9)," +
                "  treatment = gl(3,3))"
        );
        rEngine.eval("" +
                "data <- createCyclopsData(counts ~ outcome + treatment, data = dobson," +
                "  modelType = \"pr\")"
        );
        rEngine.eval("" +
                "fit <- fitCyclopsModel(data," +
                "  prior = createPrior(\"none\")," +
                "  control = createControl(noiseLevel = \"silent\"))"
        );
        rEngine.eval("" +
                "instance <- cacheCyclopsModelForJava(fit)");
        rEngine.eval("" +
                "libraryFileName <- system.file(\"libs\", \"Cyclops.so\", package = \"Cyclops\")");

        double[] mode = rEngine.eval("coef(fit)").asDoubleArray();
        int instance = rEngine.eval("instance").asInt();
        String libraryFileName = rEngine.eval("libraryFileName").asString();

        System.err.println(new WrappedVector.Raw(mode));
        System.err.println(instance);
        System.err.println(libraryFileName);

        Parameter beta = new Parameter.Default(mode.length);

        CyclopsRegressionModel model = new CyclopsRegressionModel("name",
                libraryFileName, instance, beta, true);

        System.err.println(new WrappedVector.Raw(beta.getParameterValues()));
        System.err.println(model.getLogLikelihood());

        beta.setParameterValue(0, -1.0);
        System.err.println(new WrappedVector.Raw(beta.getParameterValues()));
        System.err.println(model.getLogLikelihood());

        model.findMode();
        System.err.println(new WrappedVector.Raw(beta.getParameterValues()));
        System.err.println(model.getLogLikelihood());

    }
}
