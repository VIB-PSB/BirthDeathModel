package be.ugent.psb.setas.bdmodel.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Stack;
import java.util.TreeMap;

import be.ugent.psb.setas.bdmodel.model.ProbCalculator;
import be.ugent.psb.setas.bdmodel.parsers.NewickParser;
import be.ugent.psb.setas.bdmodel.parsers.ReadGFcountsFile;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;
import be.ugent.psb.setas.bdmodel.parsers.WGMparser;

public class LikeLihood {

	private double lambda;
	private int numOfTestRootSizes;
	private ProbCalculator probCalc;
	
	public LikeLihood(double lambda, int numOfTestRootSizes) {
		this.lambda = lambda;
		this.numOfTestRootSizes = numOfTestRootSizes;
		this.probCalc = new ProbCalculator();
	}

	public LikeLihood(double lambda, int numOfTestRootSizes,
			ProbCalculator probCalc) {
		this.lambda = lambda;
		this.numOfTestRootSizes = numOfTestRootSizes;
		this.probCalc = probCalc;
	}

	/**
	 * @param args
	 *            arln = arraylist of nodes in the post order : left-right-root
	 */


	public double getNumOfTestRootSizes() {
		return numOfTestRootSizes;
	}

	public void setNumOfTestRootSizes(int numOfTestRootSizes) {
		this.numOfTestRootSizes = numOfTestRootSizes;
	}

	public double getLambda() {
		return lambda;
	}

	public void setLambda(double lambda) {
		this.lambda = lambda;
	}

	public double lk_leaf(Node leaf, int slotInLkArray) {

		double lk = this.probCalc.probCalc(lambda, leaf.getbLen(),
				slotInLkArray, leaf.getValue());

		if (lk < 0) {
			System.err.println("leaf-node-name:   " + leaf.getName() + "\t branch-lenght:   "
					+ leaf.getbLen() + "\t value:   " + leaf.getValue() + "   lambda   "
					+ lambda);
		}
		return lk;

	}

	public double lk_normal(Node normalNode, int slotInLkArray) {

		double lk = 0;

		for (int k = 1; k < numOfTestRootSizes; k++) {

			lk += this.probCalc.probCalc(lambda, normalNode.getbLen(),
					slotInLkArray, k) * normalNode.getLkValue(k);
		}

		if (lk < 0) {
			System.err.println("normal node:   " + normalNode.getName()
					+ "   t   " + normalNode.getbLen() + "   lambda   "
					+ lambda);
		}

		return lk;
	}

	public double lk_wgm(Node wgmNode, int slotInLkArray) { // calculating lk for a WGD node as a child of another node

		double lk = 0;
		int factor = wgmNode.multFactor;
//		double factor = wgmNode.multFactor;

		for (int k = factor; k < numOfTestRootSizes; k = k + factor) {
//		for (int k = 1; k < numOfTestRootSizes; k++) {
		
			lk += this.probCalc.probCalc(lambda, wgmNode.getbLen(),
					slotInLkArray, k / factor) * wgmNode.getLkValue(k);
			
//			int floorSizeBeforeWGD =(int)(Math.floor(k/factor));	
//			if(floorSizeBeforeWGD > numOfTestRootSizes){
//				lk =0;
//			}		
//			else{
//			double remainder = (k/factor) - floorSizeBeforeWGD;
//			if(remainder >= 0.5){
//			
//			lk += this.probCalc.probCalc(lambda, wgmNode.getbLen(),
//					slotInLkArray, (int)(Math.ceil(k / factor))) * wgmNode.getLkValue(k);
//			}
//			else{
//				lk += this.probCalc.probCalc(lambda, wgmNode.getbLen(),
//						slotInLkArray, floorSizeBeforeWGD) * wgmNode.getLkValue(k);
//			}
//			}
		}

		if (lk < 0) {
			System.err.println("wgm node:   " + wgmNode.getName() + "\t"
					+ wgmNode.getbLen() + "   lambda   " + lambda);
		}

		return lk;

	}

	public double lk_vn(Node virtualNode, int slotInLkArray) {
		double lk = 0;

		for (int k = 1; k < numOfTestRootSizes; k++) {
			lk += this.probCalc.probCalc(lambda, virtualNode.getbLen(),
					slotInLkArray, k) * virtualNode.getLkValue(k);
		}

		if (lk < 0) {
			System.err.println("virtual node   " + virtualNode.getName() + "\t"
					+ virtualNode.getbLen() + "   lambda   " + lambda);
		}
		return lk;
	}

	public double decideForOneNode(Node n, int slotInLkArray) {

		if (n.isLeaf) {

			return (lk_leaf(n, slotInLkArray));
		}

		else if (n.isWGM) {

			return (lk_wgm(n, slotInLkArray));
		}

		else if (n.isVirtualNode) {
			return (lk_vn(n, slotInLkArray));
		}

		else {
			return (lk_normal(n, slotInLkArray));
		}
	}

	/* Only for a tree with 3 or more nodes */
	public double[] calcInternalLk(List<Node> arln) {

		Stack<Node> st = new Stack<Node>();
		st.add(arln.get(0));

		for (Node n : arln) {

			if (n.depth >= st.peek().depth) {
				st.add(n);
			}

			else {
				/*
				 * CASE 1: WHEN GIVEN NODE IS NOT WGM , (So it is normal, for
				 * leaves the depth is higher and it goes to the previous "if")
				 */
				if ((!n.isWGM) && (!n.isVirtualNode)) {

					Node right = st.pop();
					Node left = st.pop();

					for (int j = 1; j < numOfTestRootSizes; j++) {
						double l = decideForOneNode(left, j);
						double r = decideForOneNode(right, j);
						n.setLkValue(j, r * l);
					}
					st.add(n);
				}

				/* CASE 2: WHEN GIVEN NODE IS WGM ITSELF __>IN ROLE OF PARENT */
				else if (n.isWGM) {

					int factor = n.multFactor;

					Node child = st.pop();

					for (int j = factor; j < numOfTestRootSizes; j = j + factor) {
						double lk = decideForOneNode(child, j);
						n.setLkValue(j, lk);
					}
					st.add(n);
				}
				/*
				 * CASE 3: WHEN GIVEN NODE IS VirtualNode ITSELF __>IN ROLE OF
				 * PARENT
				 */
				else {
					Node child = st.pop();
					for (int j = 1; j < numOfTestRootSizes; j++) {
						double lk = decideForOneNode(child, j);
						n.setLkValue(j, lk);
					}
					st.add(n);
				}
			}
		}

		double[] likelihoodRoot = st.peek().getLkArray();

//		 double[] convertedLk = IncludeAngioSperms.convertLkArray(
//		 likelihoodRoot, lambda, 0.260, 0.275776, this.probCalc);
//		 double[] convertedLk =
//		 IncludeAngioSperms.convertLkArray(likelihoodRoot, lambda, 0.155383,
//		 0.275776, this.probCalc, 3);
//		 double[] convertedLk = IncludeAngioSperms.convertLkArray(
//				 likelihoodRoot, lambda, 0.222, 0.323803, this.probCalc, 2);
		return likelihoodRoot;
//		 return convertedLk;
	}

	public double[] calcInternalLk_combLogLk(Node root,
			List<List<Integer>> gfCounts, int startGF, int endGF,
			List<Node> arln) {

		double[] lks = new double[numOfTestRootSizes];
		double[] sumLoglks = new double[numOfTestRootSizes];

		for (int i = startGF; i <= endGF; i++) {

			/**
			 * These lines have to be here inside the loop, every time the GF
			 * counts change, the aryln should change, since leaves are having
			 * different values.
			 */
			SpeciesTreeParser.setLeavesValues(root, gfCounts, i);
			arln = SpeciesTreeParser.setMaxNodeSize(root, 100);

			/** for combinedLoglk All MGF we have to calculate the logarithm
			 * here, not afterwards in f in CPM, because we want the maximum of
			 * sum(loglks) not the log(sum likelihoods), in other words we can
			 * just sum up because we assume we are dealing with log values, but
			 * we can not sum up likelihoods! previously was OK because we were
			 * dealing with only 1 GF
			 **/
			lks = calcInternalLk(arln);
			double[] logLks = MathOperations.giveLogArray(lks);

			for (int j = 1; j < sumLoglks.length; j++) {
				sumLoglks[j] += logLks[j];
			}
		}
		return sumLoglks;
	}

	// Since we change lks at the end of previous method to include angiosperms,
	// this method should be used very cautiously ! maybe better not!

	// public double calcLkForRootSize(Node root, int rootSize) {
	//
	// Node left = root.getLeftChild();
	// Node right = root.getRightChild();
	//
	// ArrayList<Node> arLeft = left.postOrder(left);
	// ArrayList<Node> arRight = right.postOrder(right);
	//
	// /* To set the Likelihood array for these nodes */
	// calcInternalLk(arRight);
	// calcInternalLk(arLeft);
	//
	// double leftLk = decideForOneNode(left, rootSize);
	// double rightLk = decideForOneNode(right, rootSize);
	//
	// double lk = leftLk * rightLk;
	//
	// if (lk < 0 || lk > 1) {
	// System.err.println("Likelihood is not valid: " + lk);
	// lk = 0;
	// }
	//
	// return lk;
	// }

	// public static void main(String args[]) {
	//
	// Node root = SpeciesTreeParser.buildAndPartitionTree(args[0],args[1], 0.1,
	// 100);
	//
	// ReadGFcountsFile rgf = new ReadGFcountsFile();
	// List<List<Integer>> nonRepetetiveCounts = rgf.read_unique(args[2]);
	// ArrayList<String> gfIDs = rgf.getGfIDs_unique();
	//
	// // double wgtBlen = 0.260;
	// // double oldrootBlen = 0.275776;
	// // int numberOfObservations =1000;
	//
	// int gf = Integer.parseInt(args[3]);
	//
	// System.out.println(gfIDs.get(gf));
	//
	// SpeciesTreeParser.setLeavsValues(root, nonRepetetiveCounts, gf);
	// ArrayList<Node> speciesTree =
	// SpeciesTreeParser.setMaxNodeSizeAccordingToGF(root, nonRepetetiveCounts,
	// gf);
	//
	// for (int rootSize = 1; rootSize <= 20; rootSize++) {
	//
	// ProbCalculator probCache = new ProbCalculator();
	//
	// double[] likelihoods = new double[100];
	// double[] logLks = new double[100];
	//
	// double maxLambda = 10;
	//
	// double stepSize = 0.1;
	// int numOfLambdas = (int)(maxLambda / stepSize);
	//
	// double[] logLksOfLambdas = new double[numOfLambdas];
	//
	//
	// for (int i = 1; i < numOfLambdas; i++) {
	//
	// double lambda = i * (stepSize);
	// LikeLihood lk = new LikeLihood(lambda, 100, probCache);
	//
	// likelihoods = lk.calcInternalLk(speciesTree);
	// logLks = Log.giveLogArray(likelihoods);
	//
	// logLksOfLambdas[i] = logLks[rootSize];
	// // System.out.println(lambda+"\t"+logLksOfLambdas[i]);
	// }
	//
	// double maxLogLk = FindMaxArray.findMaxValue(logLksOfLambdas);
	// int indexMaxLk = FindMaxArray.findIndexOfMaxValue(logLksOfLambdas);
	//
	// double lambdaOpt = indexMaxLk * (stepSize);
	//
	// // assume for our CPM be.ugent.psb.setas.bdmodel.test we fix the root size to 1
	// System.out.println(lambdaOpt + "\t" + maxLogLk);
	//
	// }
	// }

	// int[] tpp = new int[root.getNumberOfLeaves()];
	// tpp[0] = 7;
	// tpp[1] = 22;
	// tpp[2] = 12;
	// tpp[3] = 7;
	// tpp[4] = 10;
	// tpp[5] = 7;
	// tpp[6] = 11;
	// tpp[7] = 5;
	// tpp[8] = 5;
	// tpp[9] = 9;
	// int[] idealObs = new int[root.getNumberOfLeaves()];
	// idealObs[0] = 10;
	// idealObs[1] = 20;
	// idealObs[2] = 10;
	// idealObs[3] = 6;
	// idealObs[4] = 8;
	// idealObs[5] = 5;
	// idealObs[6] = 10;
	// idealObs[7] = 5;
	// idealObs[8] = 5;
	// idealObs[9] = 5;

	// int numTestSizes = 7;
	// double[] lambdas = new double[16];
	// double[] lns = new double[lambdas.length];

	// lambdas[0] = 0.1;
	// lambdas[1] = 0.3;
	// lambdas[2] = 0.5;
	// lambdas[3] = 0.7;
	// lambdas[4] = 0.9;
	// lambdas[5] = 1;
	// lambdas[6] = 1.1;
	// lambdas[7] = 1.2;
	// lambdas[8] = 1.3;
	// lambdas[9] = 1.4;
	// lambdas[10] = 1.5;
	// lambdas[11] = 1.6;
	// lambdas[12] = 1.7;
	// lambdas[13] = 1.8;
	// lambdas[14] = 1.9;
	// lambdas[15] = 2;

	// lns[3] corresponds to lambda =1
	// System.out.println("for lambda=1 lnLk is: "+ lns[3]);

	// double maxLnLk = fma.ValueFindMax(lns);
	// int maxIndex = fma.IndiceFindMax(lns);
	// double bestLambda = lambdas[maxIndex];

	// }
}
