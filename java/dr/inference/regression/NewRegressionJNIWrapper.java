/*
 * NewRegressionJNIWrapper.java
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

/**
 * @author Marc A. Suchard
 */
public class NewRegressionJNIWrapper {

    /**
     * private constructor to enforce singleton instance
     */
    private NewRegressionJNIWrapper() { }

    public native double getLogLikelihood(int instance);
    
    public native void getLogLikelihoodGradient(int instance, double[] gradient);

    public native double getLogPrior(int instance);

    public native double getBeta(int instance, int index);

    public native void getBeta(int instance, double[] beta);

    public native int getBetaSize(int instance);

    public native double getHessian(int instance, int index1, int index2);

    public native void setBeta(int instance, int index, double value);

    public native void setBeta(int instance, double[] values);

    public native double getHyperprior(int instance, int index);

    public native void setHyperprior(int instance, int index, double value);

    public native void findMode(int instance);

    public native int getUpdateCount(int instance);

    public native int getLikelihoodCount(int instance);

    public native void setPriorType(int instance, int type);

    public native void makeDirty(int instance);

    public static NewRegressionJNIWrapper getCyclops(String libraryFileName)
            throws UnsatisfiedLinkError {

        if (INSTANCE == null) {
            System.err.println("Trying to load the Cyclops library...");
            System.load(libraryFileName);
            INSTANCE = new NewRegressionJNIWrapper();
            System.err.println("Cyclops library loaded.");
        }

        return INSTANCE;
    }

    private static NewRegressionJNIWrapper INSTANCE = null;
}
