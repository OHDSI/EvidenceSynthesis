package org.ohdsi.simpleDesign;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
abstract class SimpleLinearModelGradient implements GradientWrtParameterProvider {

    final SimpleLinearModel linearModel;
    final Parameter argument;
    final Parameter effect;
    final Parameter precision;
    final DesignMatrix designMatrix;

    public SimpleLinearModelGradient(SimpleLinearModel linearModel) {
        this.linearModel = linearModel;
        this.argument = linearModel.getArgument();
        this.effect = linearModel.getEffects();
        this.precision = linearModel.getPrecision();
        this.designMatrix = linearModel.getDesignMatrix();
    }

    @Override
    public Likelihood getLikelihood() {
        return linearModel;
    }

    double[] computeGradient() {

        double[] gradient = new double[argument.getDimension()];
        double[] innerProduct = linearModel.getInnerProduct();

        double tau = precision.getParameterValue(0);  // TODO Can generalize

        for (int i = 0; i < gradient.length; ++i) {
            gradient[i] = (innerProduct[i] - argument.getParameterValue(i)) * tau;
        }
        return gradient;
    }

    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this,
                getParameter().getBounds().getLowerLimit(0), getParameter().getBounds().getUpperLimit(0),
                1E-3);
    }
}
