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
package org.ohdsi.simpleDesign;

import dr.inference.model.*;

public class GroupedLinearModel extends SimpleLinearModel {

    public GroupedLinearModel(String name,
                              Parameter argument,
                              DesignMatrix designMatrix,
                              Parameter effects,
                              Parameter precision) {

        super(name, argument, designMatrix, effects, precision);
    }

    @Override
    protected double calculateLogLikelihood() {

        checkInnerProduct();

//        double tau = precision.getParameterValue(0);
//        double sd = 1.0 / Math.sqrt(tau);
//
        double logLikelihood = 0.0;
//        for (int i = 0; i < argument.getDimension(); ++i) {
//            logLikelihood += NormalDistribution.logPdf(argument.getParameterValue(i), innerProduct[i], sd);
//        }

        return logLikelihood;
    }

}
