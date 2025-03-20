package utils.bdmodel;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Stack;

public class Node {

	public Node() {
		// TODO Auto-generated constructor stub
	}

	private Node leftChild;
	private Node rightChild;
	private Node parent;

	/**
	 * If a node has no children, it is called a leaf and the boolean isleaf is True
	 * for it.
	 */
	public boolean isRoot;
	public boolean isLeaf;

	/**
	 * If this node is not a real node in the phylogeny, but a Whole Genome
	 * Duplication, we assume it is a new node and add it to the tree. The boolean
	 * isWGD is true for this kind of nodes.
	 */
	public boolean isWGD;
	public boolean isWGT;
	public boolean isWGM;
	public boolean isVirtualNode;

	public int multFactor;
	// public double multFactor;

	private String name;

	/** the branch length between a node and its parent */
	private double bLen;
	private double disToRoot;

	private int value;

	public int getValue() {
		return this.value;
	}

	public void setValue(int value) {
		this.value = value;
	}

	public int depth;

	private int maxNodeSize;
	private double[] lk = new double[maxNodeSize + 1];
	public double Lksum;

	private int optimalSize;

	private ArrayList<Node> leaves = new ArrayList<Node>();

	/**
	 * The array reference likelihood saves likelihoods of the value of a node being
	 * 1 to 100, given an array of observations from the input file (= reference
	 * observation).
	 */

	private NodeType type;

	public enum NodeType {
		INTERNAL, LEAF, WGD, VIRTUAL;
	}

	public double[] getLkArray() {
		return this.lk;
	}

	public double getLkValue(int nodeSize) {
		return this.lk[nodeSize];
	}

	public void setLkValue(int nodeSize, double likelihood) {
		if (nodeSize > this.lk.length) {
			nodeSize = this.maxNodeSize;
		}
		lk[nodeSize] = likelihood;
	}

	public int getmaxNodeSize() {
		return maxNodeSize;
	}

	public void setmaxNodeSize(int m) {
		this.maxNodeSize = m;
		double[] newLk = new double[this.maxNodeSize + 1];
		for (int i = 0; i < lk.length; i++) {
			newLk[i] = lk[i];
		}
		this.lk = newLk;
	}

	public int getOptimalSize() {
		return optimalSize;
	}

	public void setoptimalSize(int m) {
		this.optimalSize = m;
	}

	public Node getLeftChild() {
		return leftChild;
	}

	public void setLeftChild(Node leftChild) {
		this.leftChild = leftChild;
	}

	public Node getRightChild() {
		return rightChild;
	}

	public void setRightChild(Node rightChild) {
		this.rightChild = rightChild;
	}

	public Node getParent() {
		return parent;
	}

	public void setParent(Node parent) {
		this.parent = parent;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public double getbLen() {
		return bLen;
	}

	/**
	 * @return the branch length between a node and its parent
	 */
	public void setbLen(double blen) {
		this.bLen = blen;
	}

	public void setDistToRoot(double disToRoot) {

		this.disToRoot = disToRoot;
	}

	@Override
	public String toString() {
		return "Node: " + this.getName();
	}

	public int getNumberOfLeaves() {
		return leaves.size();
	}

	/** Outputs Post-order: left-right-(root) */
	public ArrayList<Node> getLeaves() {
		Stack<Node> stack = new Stack<Node>();
		Stack<Node> leaveSt = new Stack<Node>();

		if (this != null) {
			stack.push(this);

			while (!stack.empty()) {
				Node root = stack.pop();

				if (root.getLeftChild() == null && root.getRightChild() == null) {
					leaveSt.push(root);
				}

				if (root.getLeftChild() != null) {
					stack.push(root.getLeftChild());
				}
				if (root.getRightChild() != null) {
					stack.push(root.getRightChild());
				}
			}
		}

		// Changes the order reverse again to post order: desired
		leaves = new ArrayList<Node>();
		while (!leaveSt.empty()) {
			this.leaves.add(leaveSt.pop());
		}
		return this.leaves;
	}

	/**
	 * The species in GF file should be in the same order as in newick format file,
	 * i.e. post order
	 */
	public void setLeafValues(int[] observation) {
		if (leaves.isEmpty()) {
			this.getLeaves();
		}

		if (observation.length == this.getNumberOfLeaves()) {
			for (int i = 0; i < observation.length; i++) {
				this.leaves.get(i).setValue(observation[i]);

			}
		} else {
			System.out.println("Error: GF-Counts' length mismatches the number of leaves!");
		}
	}


	public void addDepthSubTree(Node n, int i) {

		if (n != null) {
			n.depth += i;
			addDepthSubTree(n.getLeftChild(), i);
			addDepthSubTree(n.getRightChild(), i);
		}
	}

	/**
	 * 
	 * @return all nodes in the tree with this node as root, traversed in postorder
	 */
	public Queue<Node> postOrder() {
		Queue<Node> q = new LinkedList<Node>();
		return postOrderRec(this, q);
	}

	private Queue<Node> postOrderRec(Node root, Queue<Node> q) {
		if (root != null) {
			postOrderRec(root.getLeftChild(), q);
			postOrderRec(root.getRightChild(), q);
			q.add(root);
		}
		return q;
	}

	/**
	 * If there are subsequent WGMs, in order for this function to work correctly,
	 * they should be mentioned in the data file in order from oldest to the
	 * youngest.
	 */
	// root-left-right
	private void addWGMRec(Node nodeInTree, String parentName, String childName, Node wgd) {

		if (nodeInTree != null) {

			if (!nodeInTree.isLeaf) {

				if (nodeInTree.getName().equals(parentName)) {

					if (nodeInTree.getLeftChild() != null && nodeInTree.getLeftChild().getName().equals(childName)) {

						addDepthSubTree(nodeInTree.getLeftChild(), 1);
						wgd.depth = (nodeInTree.depth) + 1;

						double oldBlen = nodeInTree.getLeftChild().getbLen();

						nodeInTree.getLeftChild().setParent(wgd);

						wgd.setLeftChild(nodeInTree.getLeftChild());
						wgd.getLeftChild().setbLen(oldBlen - wgd.getbLen());

						wgd.setParent(nodeInTree);
						nodeInTree.setLeftChild(wgd);
					}

					if (nodeInTree.getRightChild() != null && nodeInTree.getRightChild().getName().equals(childName)) {

						addDepthSubTree(nodeInTree.getRightChild(), 1);
						wgd.depth = (nodeInTree.depth) + 1;

						double oldBlen = nodeInTree.getRightChild().getbLen();

						nodeInTree.getRightChild().setParent(wgd);
						wgd.setRightChild(nodeInTree.getRightChild());
						wgd.getRightChild().setbLen(oldBlen - wgd.getbLen());

						wgd.setParent(nodeInTree);
						nodeInTree.setRightChild(wgd);
					}
				}

				else {
					Node left = nodeInTree.getLeftChild();
					Node right = nodeInTree.getRightChild();
					if (left != null) {
						addWGMRec(left, parentName, childName, wgd);
					}
					if (right != null) {
						addWGMRec(right, parentName, childName, wgd);
					}
				}
			}
		}

	}

	/**
	 * To read WGDs from a list and build the corresponding nodes and use addWGMRec
	 * method to add them recursively to the tree
	 * 
	 */
	public void insertWGMsToTheTree(List<List<String>> wgdList){
		/* for each line of the list = each WGD: */
		for (int i = 0; i < wgdList.size(); i++) {

			String parentName = wgdList.get(i).get(1);
			String childName = wgdList.get(i).get(2);

			Node wgm = new Node();
			wgm.isWGM = true;
			wgm.setbLen(Double.parseDouble(wgdList.get(i).get(3)));
			// Because we call it with root.insertWGM.. , so the maxNodeSize of root is
			// applied, which is set in newickparser while building tree
			wgm.setmaxNodeSize(this.getmaxNodeSize());

			if (wgdList.get(i).get(0).equalsIgnoreCase("WGD")) {
				wgm.setName("WGD-" + parentName + "-" + childName);
				wgm.isWGD = true;
				wgm.multFactor = 2;}

			if (wgdList.get(i).get(0).equalsIgnoreCase("WGT")) {
				wgm.setName("WGT-" + parentName + "-" + childName);
				wgm.isWGT = true;
				wgm.multFactor = 3;}
			addWGMRec(this, parentName, childName, wgm);
		}
	}
        
	public void insertWGMsToTheTree_Negatives(List<List<String>> wgdList){
		/* for each of the lines until keyword negative = true WGDs: */
		int trueWGDcount=0;
		for (int i = 0; i < wgdList.size(); i++) {
			if (wgdList.get(i).get(0).equalsIgnoreCase("negatives")) {
				trueWGDcount=i;
			}
		}
		// Include true positive WGDs in tree
		for (int i = 0; i < trueWGDcount; i++) {
			String parentName = wgdList.get(i).get(1);
			String childName = wgdList.get(i).get(2);

			Node wgm = new Node();
			wgm.isWGM = true;
			wgm.setbLen(Double.parseDouble(wgdList.get(i).get(3)));
			// Because we call it with root.insertWGM.. , so the maxNodeSize of root is
			// applied, which is set in newickparser while building tree
			wgm.setmaxNodeSize(this.getmaxNodeSize());

			if (wgdList.get(i).get(0).equalsIgnoreCase("WGD")) {
				wgm.setName("WGD-" + parentName + "-" + childName);
				wgm.isWGD = true;
				wgm.multFactor = 2;
			}

			if (wgdList.get(i).get(0).equalsIgnoreCase("WGT")) {
				wgm.setName("WGT-" + parentName + "-" + childName);
				wgm.isWGT = true;
				wgm.multFactor = 3;
			}
			
			addWGMRec(this, parentName, childName, wgm);
		}
	}

	public void insertWGMsToTheTree_ReverseEng(List<List<String>> wgdList, int ignoreLine){
		/* for each line of the list = each WGD: */
		for (int i = 0; i < wgdList.size(); i++) {
			if (i != ignoreLine) {
				String parentName = wgdList.get(i).get(1);
				String childName = wgdList.get(i).get(2);
				double branchLength = Double.parseDouble(wgdList.get(i).get(3));
                                
				if((i == ignoreLine+1) && (parentName.startsWith("WGD-WGD") || parentName.startsWith("WGT-WGT")) || (parentName.startsWith("WGD-WGT") || parentName.startsWith("WGT-WGD"))){
					String[] pieces = parentName.split("-");
					parentName = pieces[1].concat("-").concat(pieces[2]).concat("-").concat(pieces[3]);
					branchLength = branchLength + Double.parseDouble(wgdList.get(i-1).get(3));
				}
				else if((i == ignoreLine+1) && (parentName.startsWith("WGD") || parentName.startsWith("WGT"))){
					String[] pieces = parentName.split("-");
					parentName = pieces[1];
					branchLength = branchLength + Double.parseDouble(wgdList.get(i-1).get(3));
				}
				if((i == ignoreLine+2) && (parentName.startsWith("WGD-WGD") || parentName.startsWith("WGT-WGT")) || (parentName.startsWith("WGD-WGT") || parentName.startsWith("WGT-WGD"))){
					String[] pieces = parentName.split("-");
					parentName = pieces[0].concat("-").concat(pieces[2]).concat("-").concat(pieces[4]);
				}

				Node wg = new Node();
				wg.isWGM = true;
				wg.setbLen(branchLength);
				wg.setmaxNodeSize(this.getmaxNodeSize());

				if (wgdList.get(i).get(0).equalsIgnoreCase("WGD")) {
					wg.setName("WGD-" + parentName + "-" + childName);
					wg.isWGD = true;
					wg.multFactor = 2;
				}

				if (wgdList.get(i).get(0).equalsIgnoreCase("WGT")) {
					wg.setName("WGT-" + parentName + "-" + childName);
					wg.isWGT = true;
					wg.multFactor = 3;
				}
				addWGMRec(this, parentName, childName, wg);
			}
		}
	}
        
	public void insertWGMsToTheTree_ReverseEng_Negatives(List<List<String>> wgdList, int includeLine){
		/* for each of the lines until keyword negative = true WGDs: */
		int trueWGDcount=0;
		for (int i = 0; i < wgdList.size(); i++) {
			if (wgdList.get(i).get(0).equalsIgnoreCase("negatives")) {
				trueWGDcount=i;
			}
		}
		// Include true positive WGDs in tree
		for (int i = 0; i < trueWGDcount; i++) {
			String parentName = wgdList.get(i).get(1);
			String childName = wgdList.get(i).get(2);
			
			Node wg = new Node();
			wg.isWGM = true;
			wg.setbLen(Double.parseDouble(wgdList.get(i).get(3)));
			wg.setmaxNodeSize(this.getmaxNodeSize());

			if (wgdList.get(i).get(0).equalsIgnoreCase("WGD")) {
				wg.setName("WGD-" + parentName + "-" + childName);
				wg.isWGD = true;
				wg.multFactor = 2;
            }

            if (wgdList.get(i).get(0).equalsIgnoreCase("WGT")) {
				wg.setName("WGT-" + parentName + "-" + childName);
				wg.isWGT = true;
				wg.multFactor = 3;
            }
            addWGMRec(this, parentName, childName, wg);
		}
                
                
		// Include one true negative WGD
		// !! caution, if WGD is added on branch already containing WGDs, branch lengths may no longer make sense (in particular when adding WGD before or in between others)
                
		for (int i = trueWGDcount+1; i < wgdList.size(); i++) {
			if (i == includeLine) {
				String parentName = wgdList.get(i).get(1);
				String childName = wgdList.get(i).get(2);

				Node wg = new Node();
				wg.isWGM = true;
				wg.setbLen(Double.parseDouble(wgdList.get(i).get(3)));
				wg.setmaxNodeSize(this.getmaxNodeSize());

				if (wgdList.get(i).get(0).equalsIgnoreCase("WGD")) {
					wg.setName("WGD-" + parentName + "-" + childName);
					wg.isWGD = true;
					wg.multFactor = 2;
				}

				if (wgdList.get(i).get(0).equalsIgnoreCase("WGT")) {
					wg.setName("WGT-" + parentName + "-" + childName);
					wg.isWGT = true;
					wg.multFactor = 3;
				}
				addWGMRec(this, parentName, childName, wg);
			}
		}
	}
	
	public double calcDistToRoot() {

		double distance = this.getbLen();

		Node p = this.getParent();

		while (p.isRoot == false) {

			distance += p.getbLen();

			p = p.getParent();
		}
		return distance;
	}

	public ArrayList<Branch> findAllBranches(Node root) {

		Stack<Node> stack = new Stack<Node>();
		Stack<Branch> stackOfBranches = new Stack<Branch>();

		ArrayList<Branch> arln = new ArrayList<Branch>();

		if (root != null) {
			stack.push(root);

			while (!stack.empty()) {
				root = stack.pop();

				if (root.getLeftChild() != null) {

					Branch leftBr = new Branch(root, root.getLeftChild());

					stackOfBranches.push(leftBr);

					stack.push(root.getLeftChild());
				}
				if (root.getRightChild() != null) {

					Branch rightBr = new Branch(root, root.getRightChild());

					stackOfBranches.push(rightBr);

					stack.push(root.getRightChild());
				}
			}
		}

		while (!stackOfBranches.empty()) {
			arln.add(stackOfBranches.pop());
		}
		return arln;
	}
}