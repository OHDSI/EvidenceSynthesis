/*******************************************************************************
 * Copyright 2023 Observational Health Data Sciences and Informatics
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/

package org.ohdsi.metaAnalysis;

import dr.inference.distribution.GammaStatisticsProvider;
import dr.math.distributions.Distribution;

/**
 * A public implementation of `NormalGammaPrecisionGibbsOperator.GammaParametrization`.
 * This was originally defined a non-public static class within `dr.inference.operators.NormalGammaPrecisionGibbsOperator`.
 */
public class GammaOnPrecisionProvider implements GammaStatisticsProvider {

    private final double rate;
    private final double shape;

    /**
     * Main constructor that calculates shape and rate from a distribution's mean and variance.
     * @param mean      The mean of the distribution.
     * @param variance  The variance of the distribution.
     */
    public GammaOnPrecisionProvider(double mean, double variance) {
        if (mean == 0.0) {
            this.rate = 0.0;
            this.shape = -0.5;
        } else {
            // Standard moment matching for a Gamma distribution
            this.rate = mean / variance;
            this.shape = mean * this.rate;
        }
    }

    /**
     * Convenience constructor that takes a Distribution object.
     * @param distribution The distribution to derive parameters from.
     */
    public GammaOnPrecisionProvider(Distribution distribution) {
        this(distribution.mean(), distribution.variance());
    }

    // These methods fulfill the GammaStatisticsProvider interface contract

    public double getRate() {
        return this.rate;
    }

    public double getShape() {
        return this.shape;
    }
    public double getShape(int index) {
        return this.shape;
    }
    public double getRate(int index) {
        return this.rate;
    }
}
