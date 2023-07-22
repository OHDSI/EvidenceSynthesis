package org.ohdsi.simpleDesign;

import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class SimpleLinearModelGradientWrtEffects extends SimpleLinearModelGradientWrtArgument {

    public SimpleLinearModelGradientWrtEffects(SimpleLinearModel linearModel) {
        super(linearModel);
    }

    @Override
    public Parameter getParameter() {
        return effect;
    }

    @Override
    public int getDimension() {
        return effect.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] chain = computeGradient(); // TODO Could cache

        double[] gradient = new double[effect.getDimension()];

        for (int j = 0; j < designMatrix.getColumnDimension(); ++j) {
            double g = 0.0;
            for (int i = 0; i < designMatrix.getRowDimension(); ++i) {
                g += designMatrix.getParameterValue(i, j) * chain[i];
            }
            gradient[j] = -g;
        }

        return gradient;
    }
}
