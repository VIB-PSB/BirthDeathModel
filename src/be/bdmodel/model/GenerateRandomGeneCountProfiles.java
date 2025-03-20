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

import java.util.LinkedList;
import java.util.Queue;

/**
 * using the class generateSize, traverses the whole tree, from root to the
 * leaves, and generates sizes for each node, given the size of the root. It
 * stores the sizes of the leaves in a queue and returns it. When reaching a
 * WGD/T, it passes double/triple the generated size to continue generating
 * sizes for its descendants.
 * 
 * 
 * The class traverses the tree in post order fashion: left-right-root, so the
 * leaves are visited in this order and the generatedSizes for the leaves are in
 * this order. So when using the output of this class, one should be careful to
 * set the correct size to the correct leaf (e.g. take array of leaves also in
 * postOrder fashion)
 * 
 */
public class GenerateRandomGeneCountProfiles {

	public GenerateRandomGeneCountProfiles(int countNumberOfZeroAtLeaves, int maximumNodeSize, boolean reachNumberOfZeroCounts,
			TransitionProbabilityCalculator probabilityCalculator) {

		this.countNumberOfZeroAtLeaves = countNumberOfZeroAtLeaves;
		this.maxGeneCounts = maximumNodeSize;
		this.reachedNumberOfZero = reachNumberOfZeroCounts;
		this.generateRandomGeneCounts = new GenerateRandomGeneCount(probabilityCalculator);
	}

	public GenerateRandomGeneCountProfiles(int countNumberOfZeroAtLeaves, int maxNodeSize, boolean reachNumberOfZeroCounts,
			TransitionProbabilityCalculator probabilityCalculator, int lengthOfMCMC) {

		this.countNumberOfZeroAtLeaves = countNumberOfZeroAtLeaves;
		this.maxGeneCounts = maxNodeSize;
		this.reachedNumberOfZero = reachNumberOfZeroCounts;
		this.generateRandomGeneCounts = new GenerateRandomGeneCount(probabilityCalculator, lengthOfMCMC);
	}

	public GenerateRandomGeneCountProfiles(int countNumberOfZeroAtLeaves, int maxNodeSize, boolean reachNumberOfZeroCounts) {
		this.countNumberOfZeroAtLeaves = countNumberOfZeroAtLeaves;
		this.maxGeneCounts = maxNodeSize;
		this.reachedNumberOfZero = reachNumberOfZeroCounts;
		this.generateRandomGeneCounts = new GenerateRandomGeneCount(null);
	}

	private Queue<Integer> queue = new LinkedList<Integer>();
	private int countNumberOfZeroAtLeaves;
	private int maxGeneCounts;
	private boolean reachedNumberOfZero;
	private GenerateRandomGeneCount generateRandomGeneCounts;

	public boolean isReachedNumOfZero() {
		return reachedNumberOfZero;
	}

	public void setReachNumOfZero(boolean reachNumOfZero) {
		this.reachedNumberOfZero = reachNumOfZero;
	}

	public int getCountNumberOfZeroAtLeaves() {
		return countNumberOfZeroAtLeaves;
	}

	public void setCountNumberOfZeroAtLeaves(int countNumberOfZeroAtLeaves) {
		this.countNumberOfZeroAtLeaves = countNumberOfZeroAtLeaves;
	}

	public int getMaxNodeSize() {
		return maxGeneCounts;
	}

	public int validateRandomCounts(int randomCount) {

		if (randomCount > maxGeneCounts) {

			randomCount = maxGeneCounts;
		}

		if (randomCount < 1) {

			randomCount = 1;
		}

		return randomCount;
	}

	
	/**
	 * @param leaf
	 * @param testGeneCount size of the leave node
	 * @param lambda
	 * @return
	 */
	public int generateCountsForLeafNode(Node leaf, int testGeneCount, double lambda) {

		int generatedLeafSize = this.generateRandomGeneCounts.generateSizeForleaves(maxGeneCounts, testGeneCount, lambda, leaf.getbranchLength());
		return generatedLeafSize;

	}

	public int generateCountsForNormalNode(Node normal, int testGeneCount, double lambda) {

		int generatedSize = this.generateRandomGeneCounts.generategeneCounts(maxGeneCounts, testGeneCount, lambda, normal.getbranchLength());
		return generatedSize;

	}

	public Queue<Integer> generateQueueOfGeneCountProfiles(Node n, int testGeneCount, double lambda) {

		int testGeneCount_validated = validateRandomCounts(testGeneCount); 

		if (n.getLeftChild() != null) {

			Node leftChild = n.getLeftChild();

			if (leftChild.isLeaf) {

				int generatedGeneCountAtLeaf = generateCountsForLeafNode(leftChild, testGeneCount_validated, lambda);
				queue.add(generatedGeneCountAtLeaf);
			}

			else {
				int generatedGeneCountAtLeaf = generateCountsForNormalNode(leftChild, testGeneCount_validated, lambda);

				if (leftChild.isWGM) {

					generateQueueOfGeneCountProfiles(leftChild, leftChild.multiplicationFactor * generatedGeneCountAtLeaf, lambda);
				}

				else {
					generateQueueOfGeneCountProfiles(leftChild, generatedGeneCountAtLeaf, lambda);
				}

			}

		}

		if (n.getRightChild() != null) {

			Node rightChild = n.getRightChild();

			if (rightChild.isLeaf) {
				int generatedGeneCountAtLeaf = generateCountsForLeafNode(rightChild, testGeneCount_validated, lambda);
				queue.add(generatedGeneCountAtLeaf);
			}

			else {
				int generatedRightSize = generateCountsForNormalNode(rightChild, testGeneCount_validated, lambda);
				if (rightChild.isWGM) {

					generateQueueOfGeneCountProfiles(rightChild, rightChild.multiplicationFactor * generatedRightSize, lambda);
				} else {
					generateQueueOfGeneCountProfiles(rightChild, generatedRightSize, lambda);
				}
			}

		}
		return queue;
	}

	public static int[] queueToArray(Queue<Integer> queueOfGeneCountProfile) {
		int size = queueOfGeneCountProfile.size();
		int[] b = new int[size];

		for (int k = 0; k < size; k++) {
			b[k] = queueOfGeneCountProfile.remove();
		}
		return b;
	}

	public int[] generateGeneCountProfile(Node node, int testGeneCount, double lambda) {

		return (queueToArray(generateQueueOfGeneCountProfiles(node, testGeneCount, lambda)));
	}

}
