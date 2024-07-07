package org.ohdsi.mock;

import org.ohdsi.mcmc.Analysis;
import org.ohdsi.mcmc.Runner;
import org.ohdsi.metaAnalysis.DataModel;
import org.ohdsi.metaAnalysis.MultivariableHierarchicalMetaAnalysis;
import org.ohdsi.simpleDesign.CyclopsNormalAnalysis;

import java.util.ArrayList;
import java.util.List;

public class CyclopsHierarchicalAnalysisMock {

    public static void main(String[] args) {

        int chainLength = 1100000;
        int burnIn = 100000;
        int subSampleFrequency = 1000;

        List<MockCyclops> mocks = new ArrayList<>();
        List<DataModel> dataModels = new ArrayList<>();
        for (int i = 0; i < 4; ++i) {
            MockCyclops mock = new MockCyclops();
            mocks.add(mock);
            dataModels.add(mock.getModel());
        }

        MultivariableHierarchicalMetaAnalysis.HierarchicalMetaAnalysisConfiguration cg =
                new MultivariableHierarchicalMetaAnalysis.HierarchicalMetaAnalysisConfiguration();
        cg.tauDf = 10;

        Analysis analysis = new MultivariableHierarchicalMetaAnalysis(dataModels, cg);

        Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, 666);
        runner.run();
        runner.processSamples();

        for (MockCyclops mock : mocks) {
            mock.close();
        }
    }
}
