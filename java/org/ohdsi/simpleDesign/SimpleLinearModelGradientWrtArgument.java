package org.ohdsi.simpleDesign;

import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class SimpleLinearModelGradientWrtArgument extends SimpleLinearModelGradient {

    public SimpleLinearModelGradientWrtArgument(SimpleLinearModel linearModel) {
        super(linearModel);
    }

    @Override
    public Parameter getParameter() {
        return argument;
    }

    @Override
    public int getDimension() {
        return argument.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return computeGradient();
    }
}
