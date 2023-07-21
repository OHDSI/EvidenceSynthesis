package org.ohdsi.likelihood;

import dr.inference.model.*;
import org.ohdsi.metaAnalysis.DataModel;

public class CachedModelLikelihood extends AbstractModelLikelihood {

    private final Likelihood likelihood;
    private final Parameter parameter;

    public CachedModelLikelihood(String name, DataModel empiricalDataModel) {
        super(name);
        setId(name);

        likelihood = empiricalDataModel.getLikelihood();
        parameter = empiricalDataModel.getCompoundParameter();

        addVariable(parameter);
        likelihoodKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object o, int i) {
        throw new RuntimeException("Unknown model");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int i, Variable.ChangeType changeType) {
        if (variable == parameter) {
            likelihoodKnown = false;
        } else {
            throw new RuntimeException("Unknown variable");
        }

    }

    @Override
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    @Override
    protected void acceptState() {
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = likelihood.getLogLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;

    private double logLikelihood;
    private double storedLogLikelihood;
}
