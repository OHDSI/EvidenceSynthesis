/*
 * RegressionInCyclops.java
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
public class RegressionInCyclops {

    private final NewRegressionJNIWrapper wrapper;
    private final int instance;
    private final int dim;

    public RegressionInCyclops(String libraryFileName, int instance) {
        this.wrapper = NewRegressionJNIWrapper.getCyclops(libraryFileName);
        this.instance = instance;
        this.dim = getBetaSize();
    }

    public double getLogLikelihood() {
        return wrapper.getLogLikelihood(instance);
    }

    public void getLogLikelihoodGradient(double[] gradient) {
        assert gradient.length == dim;
        wrapper.getLogLikelihoodGradient(instance, gradient);
    }

    public double[] getLogLikelihoodGradient() {
        double[] gradient = new double[dim];
        wrapper.getLogLikelihoodGradient(instance, gradient);
        return gradient;
    }

    public double getLogPrior() {
        return wrapper.getLogPrior(instance);
    }

    public double getBeta(int index) {
        return wrapper.getBeta(instance, index);
    }

    public void getBeta(double[] beta) {
        assert beta.length == dim;
        wrapper.getBeta(instance, beta);
    }

    public int getBetaSize() {
        return wrapper.getBetaSize(instance);
    }

    public double getHessian(int index1, int index2) {
        return wrapper.getHessian(instance, index1, index2);
    }

    public void setBeta(int index, double value) {
        wrapper.setBeta(instance, index, value);
    }

    public void setBeta(double[] values) {
        assert values.length == dim;
        wrapper.setBeta(instance, values);
    }

    public double getHyperprior(int index) {
        return wrapper.getHyperprior(instance, index);
    }

    public void setHyperprior(int index, double value) {
        wrapper.setHyperprior(instance, index, value);
    }

    public void findMode() {
        wrapper.findMode(instance);
    }

    public int getUpdateCount() {
        return wrapper.getUpdateCount(instance);
    }

    public int getLikelihoodCount() {
        return wrapper.getLikelihoodCount(instance);
    }

    public void setPriorType(int type) {
        wrapper.setPriorType(instance, type);
    }

    public void makeDirty() {
        wrapper.makeDirty(instance);
    }
}
