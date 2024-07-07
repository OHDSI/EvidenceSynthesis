/*
 * MockCyclopsUsingJri.java
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

package org.ohdsi.mock;

import dr.inference.model.Parameter;
import dr.inference.regression.CyclopsRegressionModel;
import dr.inference.regression.CyclopsRegressionModelGradient;
import dr.math.matrixAlgebra.WrappedVector;
import org.rosuda.JRI.Rengine;

/**
 * @author Marc A. Suchard
 */
public class MockCyclops {

    private final CyclopsRegressionModel model;

    private static Rengine singleEngine;

    private static Rengine getEngine() {
        if (singleEngine == null) {
            singleEngine = new Rengine(new String[]{"--no-save"}, false, null);
        }
        return singleEngine;
    }

    public MockCyclops() {

        Rengine rEngine = getEngine();
//        rEngine = new Rengine(new String[]{"--no-save"}, false, null);

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

        model = new CyclopsRegressionModel("name",
                libraryFileName, instance, new Parameter.Default(mode.length), true);
    }

    public CyclopsRegressionModel getModel() {
        return model;
    }

    public void close() {
        if (singleEngine != null) {
            singleEngine.end();
        }
        singleEngine = null;
    }

    public static void main(String[] args) {

        MockCyclops mock = new MockCyclops();
        CyclopsRegressionModel model = mock.getModel();
        Parameter beta = model.getParameter();

        System.err.println(new WrappedVector.Raw(beta.getParameterValues()));
        System.err.println(model.getLogLikelihood());

        beta.setParameterValue(0, -1.0);
        System.err.println(new WrappedVector.Raw(beta.getParameterValues()));
        System.err.println(model.getLogLikelihood());

        model.findMode();
        System.err.println(new WrappedVector.Raw(beta.getParameterValues()));
        System.err.println(model.getLogLikelihood());

        CyclopsRegressionModelGradient modelGradient = new CyclopsRegressionModelGradient(model, beta);
        double[] gradient = modelGradient.getGradientLogDensity();
        System.err.println(new WrappedVector.Raw(gradient));

        mock.close();
    }
}