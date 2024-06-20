package org.ohdsi.mock;

import org.ohdsi.mcmc.Analysis;
import org.ohdsi.mcmc.Runner;
import org.ohdsi.simpleDesign.CyclopsNormalAnalysis;

public class CyclopsNormalAnalysisMock {

    public static void main(String[] args) {

        int chainLength = 1100000;
        int burnIn = 100000;
        int subSampleFrequency = 1000;

        MockCyclops mock = new MockCyclops();

        Analysis analysis = new CyclopsNormalAnalysis(mock.getModel(), 0, 10);

        Runner runner = new Runner(analysis, chainLength, burnIn, subSampleFrequency, 666);
        runner.run();
        runner.processSamples();

        mock.close();
    }
}
