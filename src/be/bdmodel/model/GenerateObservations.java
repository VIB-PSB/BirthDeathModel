package be.ugent.psb.setas.bdmodel.model;

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
	 * WGM, it passes double/triple the generated size to contiune generating
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

	public boolean isReachNumOfZero() {
		return reachNumOfZero;
	}

	public void setReachNumOfZero(boolean reachNumOfZero) {
		this.reachNumOfZero = reachNumOfZero;
	}

	public int getCountNumberOfZeroAtLeaves() {
		return countNumberOfZeroAtLeaves;
	}

	public void setCountNumberOfZeroAtLeaves(int countNumberOfZeroAtLeaves) {
		this.countNumberOfZeroAtLeaves = countNumberOfZeroAtLeaves;
	}

	public int getMaxNodeSize() {
		return maxNodeSize;
	}

	public int validateTestSize(int testSize) {

		if (testSize > maxNodeSize) {

			testSize = maxNodeSize-1;}

		if (testSize < 1) {

			testSize = 1;}

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

		int testSize = validateTestSize(testParentSize); // To avoid generating 0 in the middle of the queue.

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

	
//	/* These methods were used in SC class because VirtualNodes were not yet insterted as new nodes to the tree structure*/
//	public int genSizeLeafNode_partition(Node leaf, int testSize, double lambda) {
//
//		int generatedLeafSize = GenerateSizeWithPartition
//				.generateSizeForleaves_partitioning(maxNodeSize, testSize,
//						lambda, leaf.getbLen());
//		return generatedLeafSize;
//
//	}
//
//	public int genSizeNormalNode_partition(Node leaf, int testSize,
//			double lambda) {
//
//		int generatedLeafSize = GenerateSizeWithPartition
//				.generateSize_partitioning(maxNodeSize, testSize, lambda,
//						leaf.getbLen());
//		return generatedLeafSize;
//
//	}

//	public int genSizeAfterWGMwithRetention_partition(Node WGM, int sizeBeforeWGM,
//			double qWGD, double qWGT) {
//
//		int newParentSize = sizeBeforeWGM;
//
//		if (WGM.isWGD) {
//			newParentSize = (int) (qWGD * sizeBeforeWGM);
//		}
//
//		else if (WGM.isWGT) {
//			newParentSize = (int) (qWGT * sizeBeforeWGM);
//		}
//
//		return newParentSize;
//
//	}

//	public int genNewParentSize_NodeType(Node n, int oldParentSize,
//			double qWGD, double qWGT, double lambda, boolean reachNumOfZeros) {
//
//		int newGenSize = oldParentSize;
//
//		if (n.isLeaf) {
//
//			if (!reachNumOfZeros) {
//				newGenSize = genSizeLeafNode_partition(n, oldParentSize, lambda);
//			}
//
//			else { /* Don't generate zeros anymore */
//				newGenSize = genSizeNormalNode_partition(n, oldParentSize,
//						lambda);
//			}
//		}
//
//		else if (n.isWGM) {
//
//			int genSizeBeforeWGM = genSizeNormalNode_partition(n,
//					oldParentSize, lambda);
//			newGenSize = genSizeAfterWGMwithRetention_partition(n, genSizeBeforeWGM, qWGD,
//					qWGT);
//		}
//
//		else {
//
//			newGenSize = genSizeNormalNode_partition(n, oldParentSize, lambda);
//
//		}
//
//		return newGenSize;
//	}

//	public Queue<Integer> generateQofObservation_SC(Node n, int testParentSize,
//			double lambda, int maxZeros, double qWGD, double qWGT) {
//
//		/* The difference with generateQofObservation
//		 * There is a maximum number of zeros that can be generated at the leaves,
//		 *  we use partitioning branches up to t =0.1, assuming that the tree is not partitioned yet,
//		 * i.e. the Virtual Nodes are not part of the tree structure. 
//		 *
//		 */
//
//		int testSize = validateTestSize(testParentSize);
//
//		if (n.getLeftChild() != null) {
//
//			Node leftChild = n.getLeftChild();
//
//			int generatedSize = genNewParentSize_NodeType(leftChild, testSize,
//					qWGD, qWGT, lambda, reachNumOfZero);
//
//			if (generatedSize == 0) {
//				countNumberOfZeroAtLeaves += 1;
//			}
//
//			if (countNumberOfZeroAtLeaves > maxZeros) {
//				reachNumOfZero = true;
//			}
//
//			if (leftChild.isLeaf) {
//				queue.add(generatedSize);
//			}
//
//			else {
//				generateQofObservation_SC(leftChild, generatedSize, lambda,
//						maxZeros, qWGD, qWGT);
//			}
//
//		}
//
//		if (n.getRightChild() != null) {
//
//			Node rightChild = n.getRightChild();
//
//			int generatedSize = genNewParentSize_NodeType(rightChild, testSize,
//					qWGD, qWGT, lambda, reachNumOfZero);
//
//			if (generatedSize == 0) {
//				countNumberOfZeroAtLeaves += 1;
//			}
//
//			if (countNumberOfZeroAtLeaves > maxZeros) {
//				reachNumOfZero = true;
//			}
//
//			if (rightChild.isLeaf) {
//				queue.add(generatedSize);
//			}
//
//			else {
//				generateQofObservation_SC(rightChild, generatedSize, lambda,
//						maxZeros, qWGD, qWGT);
//			}
//
//		}
//
//		return queue;
//	}

//	public int[] generateObservation_SC(Node n, int testSize, double lambda,
//			int maxZeros, double qWGD, double qWGT) {
//
//		Queue<Integer> aQ = new LinkedList<Integer>();
//		aQ = generateQofObservation_SC(n, testSize, lambda, maxZeros, qWGD,
//				qWGT);
//		return (queueToArray(aQ));
//	}

	// public static void main(String args[]) {
	//
	// NewickParser np = new NewickParser();
	// Node root =
	// np.buildTree("/home/setas/git/StochasticBD/src/Files/Trees/37spe-MGCF5.txt",100);
	//
	// WGMparser wgm = new WGMparser();
	// List<List<String>> wgdList = wgm
	// .readInputFile("/home/setas/git/StochasticBD/src/Files/WGMs/WGDs_SC_final_MGCF5.txt");
	// root.insertWGMsToTheTree(wgdList);
	//
	// ArrayList<Node> leaves = root.getLeaves();
	// root.setNumberOfLeaves(leaves.size());
	//
	// int numberOfLeaves = root.getNumberOfLeaves();
	//
	// Queue<Node> queue = root.postOrder();
	// ArrayList<Node> aryln = new ArrayList<Node>(queue);
	//
	// GenerateObservations gen = new GenerateObservations();
	// int[] observation = new int[numberOfLeaves];
	//
	// observation = gen.generateObservation_SC(root, 1, 10,5);
	// for(int i=0; i<observation.length;i++){
	// System.out.print(observation[i]+"\t");
	// }
	//
	// //each row for one leaf
	// // int [][] nomOfrepeats =new int [101][numberOfLeaves];
	//
	// // double avg =0;
	// //
	// // int more =0;
	// // int less= 0;
	// // int equal=0;
	//
	// // for (int times=0; times<1000; times++){
	// //
	// // int[] leafCounts = gen.generateObservation(root,5, 1);
	// //
	// // root.setLeafValues(leafCounts);
	// //
	// // LikeLihood lk = new LikeLihood();
	// // ArrayList<Node> aryln = root.postOrder(root);
	// //
	// // double [] lks = lk.calcInternalLk(aryln, 1.45);
	//
	// // if(lks[5] > 2.44E-11 ){
	// // more +=1;
	// // }
	// // else if(lks[5] < 2.44E-11 ){
	// // less +=1;
	// // }
	// //
	// // else{
	// // equal+=1;
	// //
	// // }
	// //
	// //
	// //// System.out.println(lks[5]);
	// //
	// // avg+= lks[5];
	// // //i moves on the leaves:
	// //// for(int i=0; i<leafCounts.length;i++){
	// ////
	// //// nomOfrepeats[leafCounts[i]][i]+=1;
	// //// }
	// // }
	// //
	// //
	// // System.out.println("more: "+more+" less: "+less+" equal: "+equal);
	// //
	// // avg = avg/1000;
	// // System.out.println(avg);
	// // int [] modsOfLeaves= new int[numberOfLeaves];
	//
	//
	// // for(int col=0; col<numberOfLeaves;col++){
	// //
	// // int maxColValue =0;
	// // int maxColIndex=0;
	// //
	// // for(int i=0;i<101;i++){
	// //
	// // if(nomOfrepeats[i][col] > maxColValue){
	// //
	// // maxColValue = nomOfrepeats[i][col];
	// // maxColIndex=i;
	// // }
	// // }
	// //
	// // modsOfLeaves[col] = maxColIndex;
	// //
	// // System.out.println(maxColIndex);
	// // }
	//
	// }

}
