package org.ohdsi.simpleDesign;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
abstract class SimpleLinearModelGradient implements GradientWrtParameterProvider {

    final SimpleLinearModel linearModel;
    final Parameter argument;
    final Parameter effect;
    final Parameter precision;
    final DesignMatrix designMatrix;
    final List<Integer> betaToTauIndexMap;

    public SimpleLinearModelGradient(SimpleLinearModel linearModel) {
        this.linearModel = linearModel;
        this.argument = linearModel.getArgument();
        this.effect = linearModel.getEffects();
        this.precision = linearModel.getPrecision();
        this.designMatrix = linearModel.getDesignMatrix();
        this.betaToTauIndexMap = linearModel.getBetaToTauIndexMap();
    }

    @Override
    public Likelihood getLikelihood() {
        return linearModel;
    }

    /**
     * Computes the gradient component (innerProduct - argument) * precision.
     * This is an intermediate term used by subclasses to compute the full gradient
     * with respect to a specific parameter.
     *
     * @return an array holding the gradient component for each data point.
     */
    double[] computeGradient() {

        double[] gradient = new double[argument.getDimension()];
        double[] innerProduct = linearModel.getInnerProduct();

        if (betaToTauIndexMap != null) {
            // Heteroscedastic case: Use the map to find the correct precision for each argument
            for (int i = 0; i < gradient.length; ++i) {
                int precisionIndex = betaToTauIndexMap.get(i);
                double tau = precision.getParameterValue(precisionIndex);
                gradient[i] = (innerProduct[i] - argument.getParameterValue(i)) * tau;
            }
        } else {
            // Homoscedastic case: Use a single precision for all arguments
            double tau = precision.getParameterValue(0);
            for (int i = 0; i < gradient.length; ++i) {
                gradient[i] = (innerProduct[i] - argument.getParameterValue(i)) * tau;
            }
        }
        return gradient;
    }

    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this,
                getParameter().getBounds().getLowerLimit(0), getParameter().getBounds().getUpperLimit(0),
                1E0);
    }
}
