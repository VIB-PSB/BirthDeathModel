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

import java.util.List;
import java.util.Stack;

/**
 * For a given gene family, this class calculates the log-likelihoods of
 * observing those gene coutns, for a given speices tree and a given root size.
 * 
 * @param args
 *            lambda = SSD or loss rate arrayListOfNodesInOrder = arraylist of
 *            nodes in the species tree in post order : left-right-root
 * @author setas
 *
 */
public class CalculateLikeLihoods {

	private double lambda;
	private int maximumGeneCountAtNodes;
	private TransitionProbabilityCalculator probabilityCalculator;

	/**
	 * Use this method when initiating new instance of class likelihood without
	 * using cache for transition probabilities. Not recommended.
	 * 
	 * @param lambda
	 * @param maximumGeneCountAtNodes
	 *            : must be defined in the main method. It was set to 100 for the
	 *            current study.
	 */
	public CalculateLikeLihoods(double lambda, int numOfTestRootSizes) {
		this.lambda = lambda;
		this.maximumGeneCountAtNodes = numOfTestRootSizes;
		this.probabilityCalculator = new TransitionProbabilityCalculator();
	}

	/**
	 * Use this method when initiating a new instacen of class CalculateLikelihoods,
	 * using cache for transition probabilities. Recommeneded method, as it will
	 * speed up the calculations.
	 * 
	 * @param lambda
	 * @param numOfTestRootSizes
	 * @param probabilityCalculator
	 */
	public CalculateLikeLihoods(double lambda, int numOfTestRootSizes,
			TransitionProbabilityCalculator probabilityCalculator) {
		this.lambda = lambda;
		this.maximumGeneCountAtNodes = numOfTestRootSizes;
		this.probabilityCalculator = probabilityCalculator;
	}

	public double getNumOfTestRootSizes() {
		return maximumGeneCountAtNodes;
	}

	public void setNumOfTestRootSizes(int numOfTestRootSizes) {
		this.maximumGeneCountAtNodes = numOfTestRootSizes;
	}

	public double getLambda() {
		return lambda;
	}

	public void setLambda(double lambda) {
		this.lambda = lambda;
	}

	/**
	 * Calculating likelihood for a leaf node as a child of another node
	 * 
	 * @param Node
	 *            leaf
	 * @param slotInLikelihoodArray
	 * @return
	 */
	public double likelihood_leafNode(Node leaf, int slotInLikelihoodArray) {

		double likelihood = this.probabilityCalculator.probabilityCalculator(lambda, leaf.getbranchLength(),
				slotInLikelihoodArray, leaf.getGeneCountAtNode());

		if (likelihood < 0) {
			System.err.println("leaf:   " + leaf.getName() + "\t" + leaf.getbranchLength() + "\t"
					+ leaf.getGeneCountAtNode() + "   lambda   " + lambda);
		}
		return likelihood;

	}

	/**
	 * Calculating likelihood for a Normal node (Not WGM, not Virtual nor leaf Node)
	 * as a child of another node
	 * 
	 * @param
	 * @param slotInLikelihoodArray
	 * @return
	 */
	public double likelihood_normalNode(Node normalNode, int slotInLikelihoodArray) {

		double likelihood = 0;

		for (int k = 1; k < maximumGeneCountAtNodes; k++) {

			likelihood += this.probabilityCalculator.probabilityCalculator(lambda, normalNode.getbranchLength(),
					slotInLikelihoodArray, k) * normalNode.getLikelihoodValue(k);
		}

		if (likelihood < 0) {
			System.err.println("normal node:   " + normalNode.getName() + "   t   " + normalNode.getbranchLength()
					+ "   lambda   " + lambda);
		}

		return likelihood;

	}

	/**
	 * Calculating likelihood for a WGD node as a child of another node
	 * 
	 * @param wGMnode
	 * @param slotInLikelihoodArray
	 * @return
	 */
	public double likelihood_WGM(Node wGMnode, int slotInLikelihoodArray) {

		double lk = 0;
		int factor = wGMnode.multiplicationFactor;

		for (int k = factor; k < maximumGeneCountAtNodes; k = k + factor) {

			lk += this.probabilityCalculator.probabilityCalculator(lambda, wGMnode.getbranchLength(),
					slotInLikelihoodArray, k / factor) * wGMnode.getLikelihoodValue(k);
		}

		if (lk < 0) {
			System.err.println(
					"wgm node:   " + wGMnode.getName() + "\t" + wGMnode.getbranchLength() + "   lambda   " + lambda);
		}

		return lk;

	}

	/**
	 * Calculating likelihood for a Virtual node as a child of another node
	 * 
	 * @param virtualNode
	 *            : inserted on branches at equi-distances specified by int
	 *            partitionSize in the main method
	 * @param slotInLikelihoodArray
	 * @return
	 */
	public double likelihood_VirtualNode(Node virtualNode, int slotInLikelihoodArray) {
		
		double likelihood = 0;
		
		for (int k = 1; k < maximumGeneCountAtNodes; k++) {
			likelihood += this.probabilityCalculator.probabilityCalculator(lambda, virtualNode.getbranchLength(),
					slotInLikelihoodArray, k) * virtualNode.getLikelihoodValue(k);
		}
		if (likelihood < 0) {
			System.err.println("virtual node   " + virtualNode.getName() + "\t" + virtualNode.getbranchLength() + "\t"
					+ virtualNode.getGeneCountAtNode() + "   lambda   " + lambda);
		}
		return likelihood;
	}

	/**
	 * Decides whether a node is a leaf, WGM or none of these (=normal node)
	 * 
	 * @param node
	 * @param slotInLkArray
	 * @return
	 */
	public double decideForOneNode(Node node, int slotInLkArray) {

		if (node.isLeaf) {

			return (likelihood_leafNode(node, slotInLkArray));
		}

		else if (node.isWGM) {

			return (likelihood_WGM(node, slotInLkArray));
		}

		else if (node.isVirtualNode) {
			return (likelihood_VirtualNode(node, slotInLkArray));
		}

		else {
			return (likelihood_normalNode(node, slotInLkArray));
		}

	}

	/**
	 * CAUTION: Only works for a tree with 3 or more nodes ArrayList of nodes in the
	 * tree must be in post order
	 */
	public double[] calculateInternalLikelihoods(List<Node> arrayListOfNodesInOrder) {

		Stack<Node> stackOfNodes = new Stack<Node>();
		stackOfNodes.add(arrayListOfNodesInOrder.get(0));

		for (Node node : arrayListOfNodesInOrder) {

			if (node.depthOfNode >= stackOfNodes.peek().depthOfNode) {
				stackOfNodes.add(node);
			}

			else {
				if ((!node.isWGM) && (!node.isVirtualNode)) {

					Node right = stackOfNodes.pop();
					Node left = stackOfNodes.pop();

					for (int j = 1; j < maximumGeneCountAtNodes; j++) {

						double leftLikelihood = decideForOneNode(left, j);
						double rightLikelihood = decideForOneNode(right, j);

						node.setLikelihoodValue(j, rightLikelihood * leftLikelihood);
					}
					stackOfNodes.add(node);
				} else if (node.isWGM) {

					int multiplicationFactor = node.multiplicationFactor;
					Node child = stackOfNodes.pop();

					for (int j = multiplicationFactor; j < maximumGeneCountAtNodes; j = j + multiplicationFactor) {

						double likelihood = decideForOneNode(child, j);

						node.setLikelihoodValue(j, likelihood);

					}

					stackOfNodes.add(node);
				} else {
					Node child = stackOfNodes.pop();

					for (int j = 1; j < maximumGeneCountAtNodes; j++) {
						
						double lk = decideForOneNode(child, j);

						node.setLikelihoodValue(j, lk);
					}
					stackOfNodes.add(node);
				}
			}
		}
		double[] likelihoodRoot = stackOfNodes.peek().getLikelihoodArray();
		return likelihoodRoot;
	}
}
