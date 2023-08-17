package org.ohdsi.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.inference.operators.hmc.MassPreconditioningOptions;

public class HmcOptions {

    private final HamiltonianMonteCarloOperator.Options hmcOptions;
    private final MassPreconditioningOptions preconditioningOptions;
    private final MassPreconditioner preconditioner;

    public HmcOptions(GradientWrtParameterProvider gradient) {
        preconditioningOptions =  new MassPreconditioningOptions.Default(0, 0,
                0, 0, false, new Parameter.Default(1E-2), new Parameter.Default(1E+2));

        hmcOptions = new HamiltonianMonteCarloOperator.Options(0.1, 5, 0.0,
                preconditioningOptions,
                100, 1E-4,
                10, 0.1, 0.8,
                HamiltonianMonteCarloOperator.InstabilityHandler.factory("reject"));

        preconditioner = MassPreconditioner.Type.NONE.factory(gradient, null, preconditioningOptions);
    }

    public HamiltonianMonteCarloOperator.Options getHmcOptions() {  return hmcOptions; }

    public MassPreconditioningOptions getPreconditioningOptions() {  return preconditioningOptions; }

    public MassPreconditioner getPreconditioner() {  return preconditioner; }
}
