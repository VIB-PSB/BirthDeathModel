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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import be.ugent.psb.setas.bdmodel.parsers.ReadGeneFamilycountsFile;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;

/**
 * This class optimizes lambda for a given species tree, with known branch
 * length, with WGMs placed on the tree and for the gene counts of a particular
 * gene family. The optimization method which is used is called Cutting Plane
 * Methods.
 * 
 * @author setas
 *
 */
public class CuttingPlaneMethod {


	private List<Node> arrayListOfNodesOfTree;
	private int rootSize;

	private final double stepSize;
	private final double deltaLocalMoves;
	private final double toleranceDerivative;
	private final double tolerancelogLikelihood;
	private final double precisionLambda;
	private final double minimumIntervalOfLambdas;
	private final double maximumIntervalOfLambdas;
	public int counter;
	private double optimalLambda;
	private double logLikelihoodOptimalLambda;
	private TransitionProbabilityCalculator probabilityCalculator;
	private HashMap<Double, Double> logLikelihoodCache;

	public double getDeltaLocalMove() {
		return deltaLocalMoves;
	}

	public boolean isEqual(double a, double b, double PRECISION) {
		return (Math.abs(a - b) < PRECISION);
	}

	public int getRootSize() {
		return this.rootSize;
	}

	public double getOptimalLambda() {
		return optimalLambda;
	}

	public double getlogLikelihoodOptimalLambda() {
		return logLikelihoodOptimalLambda;
	}

	public double getMaxInterval() {
		return this.maximumIntervalOfLambdas;
	}

	public double getMinInterval() {
		return minimumIntervalOfLambdas;
	}

	public CuttingPlaneMethod(List<Node> speciesTree, int rootSize, double stepSize, double deltaLocalMoves,
			double tolD, double tolF, double minInterval, double maxInterval, double precision,
			TransitionProbabilityCalculator probCalc) {

		this.stepSize = stepSize;
		this.deltaLocalMoves = deltaLocalMoves;
		this.rootSize = rootSize;
		this.toleranceDerivative = tolD;
		this.tolerancelogLikelihood = tolF;
		this.minimumIntervalOfLambdas = minInterval;
		this.maximumIntervalOfLambdas = maxInterval;
		this.precisionLambda = precision;
		this.probabilityCalculator = probCalc;
		this.logLikelihoodCache = new HashMap<Double, Double>();
		this.arrayListOfNodesOfTree = speciesTree;
	}

	/**
	 * Calculates the log-likelihood of a given tree, for a given root node size,
	 * and given lambda.
	 * 
	 * @param lambda
	 * @return
	 */
	public double logLikelihood(double lambda) {
		Double value = logLikelihoodCache.get(lambda);
		if (value != null) {
			return value;
		}
		CalculateLikeLihoods likelihood = new CalculateLikeLihoods(lambda,
				arrayListOfNodesOfTree.get(0).getMaxGeneCountAtNode() + 1, probabilityCalculator);
		double[] lks = likelihood.calculateInternalLikelihoods(arrayListOfNodesOfTree);

		double[] logLikelihoods = MathematicalOperations.giveLogarithm10Array(lks);
		value = logLikelihoods[rootSize];
		logLikelihoodCache.put(lambda, value);
		return logLikelihoods[rootSize];
	}

	public double forwardDifference(double lambda) {
		double forwardDiff = 0;

		double loglikelihoodLambda = logLikelihood(lambda);
		double loglikelihoodLambdaPlusDelta = logLikelihood(lambda + stepSize);

		forwardDiff = (loglikelihoodLambdaPlusDelta - loglikelihoodLambda) / stepSize;

		return forwardDiff;
	}

	public double backwardDifference(double lambda) {
		double backwardDiff = 0;

		double loglikelihoodLambda = logLikelihood(lambda);
		double logLikelihoodLambdaMinusDelta = logLikelihood(lambda - stepSize);

		backwardDiff = (loglikelihoodLambda - logLikelihoodLambdaMinusDelta) / stepSize;

		return backwardDiff;
	}

	/**
	 * 
	 * @param lambda:
	 *            value of SSD/loss rate for which we are calulating the
	 *            differences/derivatives.
	 * @param forwardBackward:
	 *            if equal to 2, we calculate forward difference, if equal to 3 , we
	 *            calculate backward difference.
	 * @return    lineChars[0] = lambda; lineChars[1]= f(lambda); lineChars[2]= f'+ (lambda) or  f'-(lambda)
	 */
	public double[] findLineCharacteristics(double lambda, int forwardBackward) {

		double[] lineCharacteristics = new double[3];

		lineCharacteristics[0] = lambda;
		lineCharacteristics[1] = logLikelihood(lambda);
		if (forwardBackward == 2) {
			lineCharacteristics[2] = forwardDifference(lambda);
		} else if (forwardBackward == 3) {
			lineCharacteristics[2] = backwardDifference(lambda);
		}

		return lineCharacteristics;
	}

	/**
	 * Find the intersection of two lines and project it on the x-axis
	 * @param line1
	 * @param line2
	 * @return
	 */
	public double findInterSectionLines(double[] line1, double[] line2) {

		double lambdaStar = ((line1[2] * line1[0]) - (line2[2] * line2[0]) + line2[1] - line1[1])
				/ (line1[2] - line2[2]);

		if (lambdaStar <= this.minimumIntervalOfLambdas) {
			return minimumIntervalOfLambdas;
		}

		else if (lambdaStar >= this.maximumIntervalOfLambdas) {
			return maximumIntervalOfLambdas;
		}

		return lambdaStar;

	}

	/**
	 * Calculates the lambda stars (the intersection of 2 derivative lines) mapped
	 * on the lambda axis.
	 * 
	 * @param lineStar
	 * @param line1
	 * @param line2
	 * @return [lambda*_1 , lambda*_2]
	 */
	public double[] calculateLambdaStar(double[] lineStar, double[] line1, double[] line2) {

		double[] lambdaStar = new double[2];
		double lambdaStar_1 = findInterSectionLines(lineStar, line1);
		double lambdaStar_2 = findInterSectionLines(lineStar, line2);

		lambdaStar[0] = lambdaStar_1;
		lambdaStar[1] = lambdaStar_2;

		return lambdaStar;
	}

	/**
	 * Find out the value of the three input lambdas maximizes the log likelihood
	 * function.
	 * 
	 * @param lambda1
	 * @param lambda2
	 * @param lambdaStar
	 * @return
	 */
	public double argMax(double lambda1, double lambda2, double lambdaStar) {

		if (Double.isNaN(lambda1) || Double.isNaN(lambda2) || Double.isNaN(lambdaStar)) {
			throw new RuntimeException("In Cutting Plane Method: arg max of three lambdas is returned as NaN");
		}

		final double loglikelihoodLambda1 = logLikelihood(lambda1);
		final double logLikelihoodLambdaStar = logLikelihood(lambdaStar);
		final double logLikelihoodLambda2 = logLikelihood(lambda2);

		if (loglikelihoodLambda1 <= logLikelihoodLambdaStar && logLikelihoodLambda2 <= logLikelihoodLambdaStar) {
			return lambdaStar; // sub-optimality si reached
		} else if (loglikelihoodLambda1 <= logLikelihoodLambdaStar && logLikelihoodLambda2 > logLikelihoodLambdaStar) {
			return lambda2;
		}

		else if (loglikelihoodLambda1 > logLikelihoodLambdaStar && logLikelihoodLambda2 <= logLikelihoodLambdaStar) {
			return lambda1;
		}

		else if (loglikelihoodLambda1 > logLikelihoodLambdaStar && logLikelihoodLambda2 > logLikelihoodLambdaStar) {
			if (loglikelihoodLambda1 >= logLikelihoodLambda2) {
				return lambda1;
			} else {
				return lambda2;
			}
		}
		System.err.println("argMax = -1 " + lambda1 + "   " + lambda2 + "   " + lambdaStar);
		return -1;
	}

	/**
	 * Find out the index of the three input lambdas maximizes the log likelihood
	 * function.
	 * 
	 * @param lambda1
	 * @param lambda2
	 * @param lambdaStar
	 * @return
	 */
	public int argMaxIndex(double lambda1, double lambda2, double lambdaStar) {

		if (Double.isNaN(lambda1) || Double.isNaN(lambda2) || Double.isNaN(lambdaStar)) {
			throw new RuntimeException("NaN is not good!");
		}

		final double logLikelihoodambda1 = logLikelihood(lambda1);
		final double logLikelihoodLambdaStar = logLikelihood(lambdaStar);
		final double logLikelihoodLambda2 = logLikelihood(lambda2);

		if (logLikelihoodambda1 <= logLikelihoodLambdaStar && logLikelihoodLambda2 <= logLikelihoodLambdaStar) {
			return 2; // sub-optimality is reached
		}

		else if (logLikelihoodambda1 <= logLikelihoodLambdaStar && logLikelihoodLambda2 > logLikelihoodLambdaStar) {
			return 1;
		} else if (logLikelihoodambda1 > logLikelihoodLambdaStar && logLikelihoodLambda2 <= logLikelihoodLambdaStar) {
			return 0;
		} else if (logLikelihoodambda1 > logLikelihoodLambdaStar && logLikelihoodLambda2 > logLikelihoodLambdaStar) {

			if (logLikelihoodambda1 >= logLikelihoodLambda2) {
				return 0;
			} else {
				return 1;
			}
		}
		System.err.println("argMaxIndex = -1 " + lambda1 + "   " + lambda2 + "   " + lambdaStar);
		return -1;
	}

	/**
	 * 
	 * @param lineStar:
	 *            line drawn from lambda_star, the intersection of previous L1 and
	 *            L2 from the begining and end of the previous interval.
	 * @param lineMin:
	 *            line drawn for lambda_min at the minimum of the current interval.
	 * @param lineMax:
	 *            line drawn for lambda_max at the maximum of the current interval.
	 * @return array of double as described in {@link lineCharacteristics}
	 */
	public double[] calculateLineSubOptimal(double[] lineStar, double[] lineMin, double[] lineMax) {

		double[] lineSubOptimal = new double[3];

		if (lineStar[0] >= lineMax[0]) {
			return (findLineCharacteristics(lineMax[0], 3));
		}

		else if (lineStar[0] <= lineMin[0]) {
			return (findLineCharacteristics(lineMin[0], 2));
		}

		else {

			double[] lambdaStar12 = calculateLambdaStar(lineStar, lineMin, lineMax);
			double lambdaTemp = argMax(lineStar[0], lambdaStar12[0], lambdaStar12[1]);

			if (isEqual(lambdaTemp, lineStar[0], precisionLambda)) {
				return (lineStar);

			}
			else {
				if (isEqual(lambdaTemp, lambdaStar12[1], precisionLambda)) {

					double[] lineMinNew = lineStar;

					double[] lineStar2 = findLineCharacteristics(lambdaStar12[1], 2);

					lineSubOptimal = calculateLineSubOptimal(lineStar2, lineMinNew, lineMax);
				}

				if (isEqual(lambdaTemp, lambdaStar12[0], precisionLambda)) {

					double[] lineMaxNew = lineStar;

					double[] lineStar1 = findLineCharacteristics(lambdaStar12[0], 2);

					lineSubOptimal = calculateLineSubOptimal(lineStar1, lineMin, lineMaxNew);
				}
			}
		}
		return lineSubOptimal;
	}

	/**
	 * 
	 * @param lambdaSubOptimal:
	 *            lambda sub-optimal taken from the characteristics of line sub-optimal
	 *            optimal produced by calcLineSubOptimal.
	 * @param fb:
	 *            forward or backward difference
	 * @param dLocalMoves:
	 *            delta local moves
	 * @return
	 */
	public double localMoves(double lambdaSubOptimal, int forwardBackward, double deltaLocalMoves) {

		if (lambdaSubOptimal <= minimumIntervalOfLambdas || lambdaSubOptimal >= maximumIntervalOfLambdas) {
			return lambdaSubOptimal;
		} else {

			double[] lineSubOptimal = findLineCharacteristics(lambdaSubOptimal, forwardBackward);

			if (isEqual(lineSubOptimal[2], 0, toleranceDerivative)) {
				return lambdaSubOptimal;
			}

			double deltaLocalMove_signed = deltaLocalMoves;
			if (lineSubOptimal[2] < 0) {
				deltaLocalMove_signed = (-1) * deltaLocalMoves;
			}

			double lambdaNew = lambdaSubOptimal + deltaLocalMove_signed;

			final double logLikelihoodLambdaNew = logLikelihood(lambdaNew);
			final double logLikelihoodLambdaSubOptimal = logLikelihood(lambdaSubOptimal);

			if (isEqual(logLikelihoodLambdaNew, logLikelihoodLambdaSubOptimal, tolerancelogLikelihood)) {
				return lambdaSubOptimal;
			} else if ((logLikelihoodLambdaNew - logLikelihoodLambdaSubOptimal) > tolerancelogLikelihood) {
				return localMoves(lambdaNew, forwardBackward, deltaLocalMoves);
			} else if ((logLikelihoodLambdaSubOptimal - logLikelihoodLambdaNew) > tolerancelogLikelihood) {
				return localMoves(lambdaSubOptimal, forwardBackward, deltaLocalMoves * 0.5);
			}

		}
		return lambdaSubOptimal;
	}

	public void findOptimalLambda() {
		double[] lineMin = findLineCharacteristics(minimumIntervalOfLambdas, 2);
		double[] lineMax = findLineCharacteristics(maximumIntervalOfLambdas, 3);

		double lambdaStar = findInterSectionLines(lineMin, lineMax);
		double[] lineStar = findLineCharacteristics(lambdaStar, 2);
		double[] lineSubOptimal = calculateLineSubOptimal(lineStar, lineMin, lineMax);

		this.optimalLambda = localMoves(lineSubOptimal[0], 2, deltaLocalMoves);
		this.logLikelihoodOptimalLambda = this.logLikelihood(optimalLambda);

	}

	/**
	 * The main method of the program
	 * 
	 * @param args
	 *            args[0]: a text file with the Newick format of the tree in it,
	 *            ending in the root node , followed by ";".
	 * 
	 *            args[1]: a text file representing WGM events. Each WGM must be
	 *            present in one line. In case successive WGMs happened on one
	 *            branch the events should be present from the oldest to the
	 *            youngest. For example WGD from node A to B with branch length x is
	 *            presented as: WGD, A, B, x where x is the distance of the WGD node
	 *            to the parent node.
	 * 
	 *            args[2]: Gene family counts. this file must contain at every line,
	 *            a gene family ID, followed by gene counts in different species in
	 *            the same order that the species appear in the Newick format tree.
	 * 
	 *            args[3]: An integer representing the number of gene family in the
	 *            gene counts file, starting from 0.
	 */
	public static void main(String[] args) {

		final double stepSize = 1e-4;// step size of calculating the numerical derivatives
		final double deltaLocalMoves = 1e-1;
		final double tolD = 1e-3;
		final double tolF = 1e-4;
		final double minInterval = 1e-2;
		final double maxInterval = 9.9999;
		final double precisionLambda = 1e-5;
		final double partitionSizeOfBranches = 0.1;
		final int maximumGeneCount = 100;
		final int numberOfGeneratedSampleGeneCounts = 1000;

		Node root = SpeciesTreeParser.buildAndPartitionTree(args[0], args[1], partitionSizeOfBranches,
				maximumGeneCount);

		ReadGeneFamilycountsFile readGeneFamilyCountsFile = new ReadGeneFamilycountsFile();
		List<List<Integer>> genefamily_Counts = readGeneFamilyCountsFile.readGeneCountsFile(args[2]);

		ArrayList<String> genefamilyIDs = readGeneFamilyCountsFile.getGfIDs();

		int geneFamilyNumber = Integer.parseInt(args[3]);
		System.out.println(genefamilyIDs.get(geneFamilyNumber));

        int rootSize = Integer.parseInt(args[4]);

			SpeciesTreeParser.setLeavesValues(root, genefamily_Counts, geneFamilyNumber);

			ArrayList<Node> speciesTree = SpeciesTreeParser.setMaximumGeneCount(root, maximumGeneCount);

			TransitionProbabilityCalculator probabilityCalculator = new TransitionProbabilityCalculator();

			CuttingPlaneMethod cuttingPlaneMethod = new CuttingPlaneMethod(speciesTree, rootSize, stepSize,
					deltaLocalMoves, tolD, tolF, minInterval, maxInterval, precisionLambda, probabilityCalculator);

			cuttingPlaneMethod.findOptimalLambda();

			System.out.print(rootSize + "\t" + cuttingPlaneMethod.getOptimalLambda() + "\t"
					+ cuttingPlaneMethod.getlogLikelihoodOptimalLambda());

			Pvalues pvalues = new Pvalues(root, speciesTree, cuttingPlaneMethod.getOptimalLambda(),
					numberOfGeneratedSampleGeneCounts, probabilityCalculator, rootSize);

			double pValue = pvalues.calculateConditionalPvalues(rootSize, rootSize,
					cuttingPlaneMethod.getlogLikelihoodOptimalLambda());

			System.out.print("\t" + pValue + "\n");

	}

}
