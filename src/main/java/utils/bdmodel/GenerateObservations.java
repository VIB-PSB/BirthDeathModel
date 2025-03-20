package utils.bdmodel;

import java.util.LinkedList;
import java.util.Queue;

public class GenerateObservations {
		
	public GenerateObservations(int countNumberOfZeroAtLeaves, int maxNodeSize, boolean reachNumOfZero, ProbCalculator probCache) {

		this.countNumberOfZeroAtLeaves = countNumberOfZeroAtLeaves;
		this.maxNodeSize = maxNodeSize;
		this.reachNumOfZero = reachNumOfZero;
		this.gs = new GenerateSize(probCache);
	}
	
	public GenerateObservations(int countNumberOfZeroAtLeaves, int minNodeSize, int maxNodeSize, boolean reachNumOfZero, ProbCalculator probCache , int lMCMC) {

		this.countNumberOfZeroAtLeaves = countNumberOfZeroAtLeaves;
		this.maxNodeSize = maxNodeSize;
		this.reachNumOfZero = reachNumOfZero;
		this.gs = new GenerateSize(probCache, lMCMC, minNodeSize, maxNodeSize);
	}
	

	public GenerateObservations(int countNumberOfZeroAtLeaves, int maxNodeSize, boolean reachNumOfZero) {
		this.countNumberOfZeroAtLeaves = countNumberOfZeroAtLeaves;
		this.maxNodeSize = maxNodeSize;
		this.reachNumOfZero = reachNumOfZero;
		this.gs = new GenerateSize();
	}
	
	/**
	 * using the class generateSize, traverses the whole tree, from root to the
	 * leaves, and generates sizes for each node, given the size of the root. It
	 * stores the sizes of the leaves in a queue and returns it. When reaching a
	 * WGM, it passes double/triple the generated size to continue generating
	 * sizes for its descendants.
	 * 
	 * 
	 * The class traverses the tree in post order fashion: left-right-root, so
	 * the leaves are visited in this order and the generatedSizes for the
	 * leaves are in this order. So when using the output of this class, one
	 * should be careful to set the correct size to the correct leaf (e.g. take
	 * array of leaves also in postOrder fashion)
	 * 
	 * @param testsize
	 *            : size of the root node
	 */

	private Queue<Integer> queue = new LinkedList<Integer>();
	private int countNumberOfZeroAtLeaves;
	private int maxNodeSize;
	private boolean reachNumOfZero;
	private GenerateSize gs;

	public int validateTestSize(int testSize) {

		if (testSize > maxNodeSize) {
            testSize = maxNodeSize;
                }
		if (testSize < 1) {
			testSize = 1;
                }
		return testSize;
	}

	public int genSizeLeafNode(Node leaf, int testSize, double lambda) {

		int generatedLeafSize = this.gs.generateSizeForleaves(
				testSize, lambda, leaf.getbLen());
		return generatedLeafSize;

	}

	public int genSizeNormalNode(Node normal, int testSize, double lambda) {

		int generatedSize = this.gs.generateSize(testSize,
				lambda, normal.getbLen());
		return generatedSize;

	}

	public Queue<Integer> generateQofObservation(Node n, int testParentSize,
			double lambda) {

        // Useful for making sure GF sizes after WGD do not exceed maxNodeSize, O's should not be produced in code below except for leaf nodes
		int testSize = validateTestSize(testParentSize); 

		if (n.getLeftChild() != null) {

			Node leftChild = n.getLeftChild();

			if (leftChild.isLeaf) {

				int generatedLeafSize = genSizeLeafNode(leftChild, testSize,
						lambda);
				queue.add(generatedLeafSize);
			}

			else { int generatedLeftSize = genSizeNormalNode(leftChild, testSize,lambda);

				if (leftChild.isWGM) {

					generateQofObservation(leftChild, leftChild.multFactor
							* generatedLeftSize, lambda);}

				else {generateQofObservation(leftChild, generatedLeftSize, lambda);}
			}
		}

		if (n.getRightChild() != null) {

			Node rightChild = n.getRightChild();

			if (rightChild.isLeaf) {
				int generatedLeafSize = genSizeLeafNode(rightChild, testSize,
						lambda);
				queue.add(generatedLeafSize);
			}

			else {
				int generatedRightSize = genSizeNormalNode(rightChild,
						testSize, lambda);
				if (rightChild.isWGM) {

					generateQofObservation(rightChild, rightChild.multFactor
							* generatedRightSize, lambda);
				} else {
					generateQofObservation(rightChild, generatedRightSize,
							lambda);
				}
			}

		}
		return queue;
	}

	public static int[] queueToArray(Queue<Integer> q) {
		int size = q.size();
		int[] b = new int[size];

		for (int k = 0; k < size; k++) {
			b[k] = q.remove();
		}
		return b;
	}

	public int[] generateObservation(Node n, int testSize, double lambda) {

		Queue<Integer> aQ = new LinkedList<Integer>();
		aQ = generateQofObservation(n, testSize, lambda);
		return (queueToArray(aQ));
	}
}
