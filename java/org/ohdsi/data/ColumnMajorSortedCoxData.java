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
package org.ohdsi.data;

/**
 * @author Marc A. Suchard
 */
public class ColumnMajorSortedCoxData extends SortedCoxData {

    public ColumnMajorSortedCoxData(int[] y, double[] x, int[] strata, int[] n, double[] weight) {
        super(y, transpose(x, y.length, x.length / y.length), strata, n, weight);
    }

    public static double[] transpose(double[] in, int nRows, int nCols) {
        double[] out = new double[in.length];

        int destination = 0;
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                out[destination++] = in[j * nRows + i];
            }
        }

        return out;
    }
}
