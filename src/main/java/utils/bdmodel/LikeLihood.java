package utils.bdmodel;

import java.util.List;
import java.util.Stack;

import utils.parsers.SpeciesTreeParser;

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
	 * arln = arraylist of nodes in the post order: left-right-root
	 */


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

	public double lk_wgm(Node wgmNode, int slotInLkArray) { // Calculating lk for a WGD node as a child of another node

		double lk = 0;
		int factor = wgmNode.multFactor;

		for (int k = factor; k < numOfTestRootSizes; k = k + factor) {
			lk += this.probCalc.probCalc(lambda, wgmNode.getbLen(),
					slotInLkArray, k / factor) * wgmNode.getLkValue(k);
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

		return likelihoodRoot;
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
}
