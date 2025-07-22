/*******************************************************************************
 * Copyright 2025 Observational Health Data Sciences and Informatics
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
package org.ohdsi.mcmc;

import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public interface Analysis {

	Likelihood getJoint();
	OperatorSchedule getSchedule();
	List<Loggable> getLoggerColumns();

	abstract class Base implements Analysis {

		protected Likelihood likelihood;
		protected Likelihood prior;
		protected Likelihood joint;
		protected OperatorSchedule schedule;
		protected Parameter beta;

		@Override
		public List<Loggable> getLoggerColumns() {

			likelihood.setId("likelihood");
			prior.setId("prior");

			List<Loggable> columns = new ArrayList<>();
			columns.add(likelihood);
			columns.add(prior);
			columns.add(beta);

			return columns;
		}

		@Override
		public Likelihood getJoint() {
			return joint;
		}

		@Override
		public OperatorSchedule getSchedule() {
			return schedule;
		}
	}

	static Parameter makeParameter(String name, int dim) {
		Parameter param = new Parameter.Default(name, dim, 0.0);

		double[] uppers = new double[dim];
		double[] lowers = new double[dim];
		Arrays.fill(uppers, Double.POSITIVE_INFINITY);
		Arrays.fill(lowers, Double.NEGATIVE_INFINITY);

		param.addBounds(new Parameter.DefaultBounds(uppers, lowers));
		return param;
	}
}
