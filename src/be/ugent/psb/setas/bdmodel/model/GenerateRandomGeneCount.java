package be.ugent.psb.setas.bdmodel.model;

/*
 * #%L
 * BirthDeathModel
 * %%
 * Copyright (C) 2017 VIB/PSB/UGent - Setareh Tasdighian
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

import java.util.Random;

import be.ugent.psb.setas.bdmodel.model.TransitionProbabilityCalculator;

/**
 * Given the size of a Parent node, generates a size for its child node,
 * sampled from a probablity distribution described in the class
 * Probcalculator, via MCMC method, with 10000 iterations.
 */
public class GenerateRandomGeneCount {

	private TransitionProbabilityCalculator probCalc;
	private int lengthOfMCMC;

	public GenerateRandomGeneCount() {
		this.probCalc = new TransitionProbabilityCalculator();
		this.lengthOfMCMC = 10000;
	}

	public GenerateRandomGeneCount(TransitionProbabilityCalculator probabilityCalculator) {
		this.probCalc = probabilityCalculator;
		this.lengthOfMCMC = 10000;

	}

	public GenerateRandomGeneCount(TransitionProbabilityCalculator probabilityCalculator, int lenghtOfMCMC) {
		this.probCalc = probabilityCalculator;
		this.lengthOfMCMC = lenghtOfMCMC;

	}

	/**
	 * Parent nodes can not take gene count =0 since the
	 * be.ugent.psb.setas.bdmodel.model does not support 0 values in the
	 * internal nodes: absorbing state of the system
	 * 
	 * @param maximumNodeSize maximum gene counts
	 * 
	 * @param parentSize
	 *            = gene counts at the parent node
	 * @param lambda
	 *            = gene duplication/loss rate
	 * @param t
	 *            = branch length
	 *            
	 * @return Generated gene count
	 */
	public int generategeneCounts(int maximumNodeSize, int parentSize, double lambda, double branchLength) {

		if (parentSize > maximumNodeSize) {
			parentSize = maximumNodeSize;
		}

		if (parentSize < 1) {
			System.err.println("Error: Internal nodes must have positive sizes. Size changed to 1");
			parentSize = 1;
		}
		Random random = new Random();

		double probabilityOfNoChange = this.probCalc.probabilityCalculator(lambda, branchLength, parentSize,
				parentSize);
		int nextGeneratedSize = parentSize;

		for (int k = 1; k < lengthOfMCMC + 1; k++) {

			int proposedChildSize = 1 + random.nextInt(maximumNodeSize);

			double fraction = this.probCalc.probabilityCalculator(lambda, branchLength, nextGeneratedSize,
					proposedChildSize) / probabilityOfNoChange;

			double alpha = Math.min(1.0, fraction);
			double u = random.nextDouble();

			if (alpha >= u) {
				nextGeneratedSize = proposedChildSize;
				probabilityOfNoChange = this.probCalc.probabilityCalculator(lambda, branchLength, nextGeneratedSize,
						nextGeneratedSize);
			}
		}

		return nextGeneratedSize;
	}

	/**
	 * Because at the leaves, the gene family can go extinct, so we can generate 0
	 * for leaf sizes : c=0
	 */
	public int generateSizeForleaves(int maxNodeSize, int parentSize, double lambda, double branchLength) {

		Random random = new Random();
		if (parentSize > maxNodeSize) {
			parentSize = maxNodeSize;
		}

		if (parentSize < 1) {
			System.err.println("Error: Internal nodes must have positive sizes. Size changed to 1");
			parentSize = 1;
		}

		double probNoChange = this.probCalc.probabilityCalculator(lambda, branchLength, parentSize, parentSize);
		int nextGeneratedSize = parentSize;

		for (int k = 1; k < lengthOfMCMC + 1; k++) {

			int proposedChildSize = random.nextInt(maxNodeSize + 1);

			double fraction = this.probCalc.probabilityCalculator(lambda, branchLength, nextGeneratedSize,
					proposedChildSize) / probNoChange;

			double alpha = Math.min(1.0, fraction);
			double u = random.nextDouble();

			if (alpha >= u) {

				nextGeneratedSize = proposedChildSize;

				if (proposedChildSize == 0) {
					return 0;
				}

				else {
					probNoChange = this.probCalc.probabilityCalculator(lambda, branchLength, nextGeneratedSize,
							nextGeneratedSize);
				}
			}
		}

		return nextGeneratedSize;

	}

}
