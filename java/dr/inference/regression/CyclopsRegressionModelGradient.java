/*
 * CyclopsRegressionModelGradient.java
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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;

/**
 * @author Marc A. Suchard
 */
public class CyclopsRegressionModelGradient implements GradientWrtParameterProvider {

    private final CyclopsRegressionModel cyclopsModel;
    private final Parameter beta;

    public CyclopsRegressionModelGradient(CyclopsRegressionModel cyclopsModel,
                                          Parameter beta) {
        this.cyclopsModel = cyclopsModel;
        this.beta = beta;
    }

    @Override
    public Likelihood getLikelihood() { return cyclopsModel; }

    @Override
    public Parameter getParameter() { return beta; }

    @Override
    public int getDimension() { return beta.getDimension(); }

    @Override
    public double[] getGradientLogDensity() { return cyclopsModel.getGradientWrtBeta(); }
}
