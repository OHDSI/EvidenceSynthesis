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

import dr.evomodel.operators.PrecisionMatrixGibbsOperator;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.MultivariateNormalGibbsOperator;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.math.distributions.WishartDistribution;
import dr.math.distributions.WishartStatistics;
import dr.math.matrixAlgebra.IllegalDimension;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.ohdsi.metaAnalysis.MultivariableHierarchicalMetaAnalysis.diagonalMatrixParameter;
import static org.ohdsi.metaAnalysis.MultivariableHierarchicalMetaAnalysis.diagonalScaleMatrix;

@SuppressWarnings("unused")
public interface MultivariatePrior {

	List<Parameter> getParameters();

	Likelihood getPrior();

	List<MCMCOperator> getOperators(double weight, AdaptationMode mode);

	abstract class Base implements MultivariatePrior {

		final List<Parameter> betas;
		final List<Parameter> parameters;
		final List<MCMCOperator> operators;
		final List<Likelihood> distributions;

		Likelihood prior;

		private List<Parameter> getAllBetas(List<DataModel> dataModels) {
			List<Parameter> betas = new ArrayList<>();
			for (DataModel dataModel : dataModels) {
				betas.add(dataModel.getCompoundParameter());
			}

			return betas;
		}

		protected Base(List<DataModel> dataModels) {
			this.betas = getAllBetas(dataModels);
			this.parameters = new ArrayList<>();
			this.operators = new ArrayList<>();
			this.distributions = new ArrayList<>();
		}

		@Override
		public List<Parameter> getParameters() {
			return parameters;
		}

		@Override
		public Likelihood getPrior() {
			if (prior == null) {
				prior = new CompoundLikelihood(distributions);
			}
			return prior;
		}
	}

	@SuppressWarnings("unused")
	class MultivariateNormal extends  Base {

		public MultivariateNormal(List<DataModel> dataModels,
								  MultivariableHierarchicalMetaAnalysis.HierarchicalMetaAnalysisConfiguration cg) {
			super(dataModels);

			final int analysisDim = betas.get(0).getDimension();

			// Build hierarchical priors and operators
			Parameter mu = new Parameter.Default("mean", analysisDim, cg.startingMu);
			mu.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, analysisDim));

			MatrixParameter tau2 = diagonalMatrixParameter("tau", analysisDim, cg.startingTau);

			MultivariateDistributionLikelihood hierarchy = new MultivariateDistributionLikelihood(
					new MultivariateNormalDistributionModel(mu, tau2));

			for (Parameter beta : betas) {
				hierarchy.addData(beta);
			}

			double[] muPriorMean = new double[analysisDim];
			Arrays.fill(muPriorMean, cg.muMean);
			double muPriorPrecision = 1 / (cg.muSd * cg.muSd);

			MultivariateDistributionLikelihood meanPrior = new MultivariateDistributionLikelihood(
					new MultivariateNormalDistribution(muPriorMean,muPriorPrecision));
			meanPrior.addData(mu);

			MultivariateDistributionLikelihood tau2Prior = new MultivariateDistributionLikelihood(
					new WishartDistribution(cg.tauDf, diagonalScaleMatrix(analysisDim, cg.tauScale)));
			tau2Prior.addData(tau2);

			distributions.add(hierarchy);
			distributions.add(meanPrior);
			distributions.add(tau2Prior);

			parameters.add(mu);
			parameters.add(tau2);

			MCMCOperator meanOperator = null;
			try {
				meanOperator = new MultivariateNormalGibbsOperator(hierarchy, meanPrior, 1.0);
			} catch (IllegalDimension e) {
				e.printStackTrace();
			}

			MCMCOperator tau2Operator = new PrecisionMatrixGibbsOperator(hierarchy,
					(WishartStatistics) tau2Prior.getDistribution(), 1.0);

			operators.add(meanOperator);
			operators.add(tau2Operator);
		}

		@Override
		public List<MCMCOperator> getOperators(double weight, AdaptationMode mode) {
			return operators;
		}
	}

	@SuppressWarnings("unused")
	class IndependentNormal extends Base {

		public IndependentNormal(List<DataModel> dataModels,
								 MultivariableHierarchicalMetaAnalysis.HierarchicalMetaAnalysisConfiguration cg) {
			super(dataModels);

			double mean = cg.muMean;
			double sd = cg.muSd;
			DistributionLikelihood likelihood = new DistributionLikelihood(new NormalDistribution(mean, sd));

			for (Parameter p : betas) {
				likelihood.addData(p);
			}

			distributions.add(likelihood);
		}

		@Override
		public List<MCMCOperator> getOperators(double weight, AdaptationMode mode) {
			return operators;
		}
	}
}
