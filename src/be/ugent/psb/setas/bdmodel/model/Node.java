package be.ugent.psb.setas.bdmodel.model;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Stack;

import be.ugent.psb.setas.bdmodel.parsers.NewickParser;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;

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
	// private int numberOfLeaves;

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

	public double getDistToRoot() {
		return disToRoot;
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

	// public void setNumberOfLeaves(int n) {
	// this.numberOfLeaves = n;
	// }

	/**
	 * old: used in Parser/BuildTree which is used for test trees in line format
	 **/
	public Node AddChildren(String a, String b, String c, double d, double e) {
		Node parent = new Node();
		parent.setName(a);

		Node leftChild = new Node();
		Node rightChild = new Node();

		leftChild.setName(b);
		leftChild.setbLen(d);

		rightChild.setName(c);
		rightChild.setbLen(e);

		parent.setLeftChild(leftChild);
		parent.setRightChild(rightChild);

		return parent;
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

		// changes the order reverse again to post order: desired
		leaves = new ArrayList<Node>();
		while (!leaveSt.empty()) {
			this.leaves.add(leaveSt.pop());
		}
		return this.leaves;
	}

	// public void setLeaves(ArrayList<Node> l) {
	//
	// this.leaves = l;
	// }

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

	/*
	 * This method is static because we are not using any infromation of the node
	 * that is calling it
	 */
	public static Node resetAllSettings(String fileName, int maxNodeSize) {

		NewickParser np = new NewickParser();
		Node root = np.buildTree(fileName, maxNodeSize);
		root.getLeaves();

		return (root);

	}

	public static Node resetAllSettings_new(String treeFile, String WGDfile, double partitionSize, int maxNodeSize) {

		return (SpeciesTreeParser.buildInsertWGDsandPartitionTree(treeFile, WGDfile, partitionSize, maxNodeSize));

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

	public Node findNodeWithName(String name) {
		Queue<Node> allNodes = this.postOrder();
		for (Node node : allNodes) {
			if (node.getName().equals(name)) {
				return node;
			}
		}
		return null;
	}

	public void addWGM(String pName, String cName, Node w) {
		addWGMRec(this, pName, cName, w);
	}

	/**
	 * If there are subsequent WGMs, in order for this function to work correctly,
	 * they should be mentioned in the data file in order from oldest to the
	 * youngest.
	 */
	// root-leaft-right
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
	public void insertWGMsToTheTree(List<List<String>> wgdList)

	{/* for each line of the list = each WGD: */
		for (int i = 0; i < wgdList.size(); i++) {

			String parentName = wgdList.get(i).get(1);
			String childName = wgdList.get(i).get(2);

			Node wgm = new Node();
			wgm.isWGM = true;
			wgm.setbLen(Double.parseDouble(wgdList.get(i).get(3)));
			// because we call it with root.insertWGM.. , so the maxNodeSize of root is
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
			
//			System.out.println(wgm.getName()+ "\t"+ wgm.getbLen());
			addWGMRec(this, parentName, childName, wgm);
		}
	}

	public void insertWGMsToTheTree_ReverseEng(List<List<String>> wgdList, int ignoreLine) {
		/* for each line of the list = each WGD: */
		for (int i = 0; i < wgdList.size(); i++) {
			if (i != ignoreLine) {
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
	
	public void deleteNode(Node n) {

		if (n.getLeftChild() != null) {
			
			n.setLeftChild(null);
		}
		if (n.getRightChild() != null) {
			
			n.setRightChild(null);
		}

		if (n.getParent().getLeftChild()!=null && n.getParent().getLeftChild().getName().equalsIgnoreCase(n.getName())) {
			
			n.getParent().setLeftChild(null);
			
			if(n.getParent().getRightChild()!=null && !n.isWGD && !n.isWGT) { //parent was a normal node
				
				n.getParent().isVirtualNode=true; //for calculations of the likelihood

			}
			
			if(n.getParent().getLeftChild()==null && n.getParent().getRightChild()==null) { //parent was virtual or WGM
				
				deleteNode(n.getParent());
			}
		}

		if (n.getParent().getRightChild()!=null && n.getParent().getRightChild().getName().equalsIgnoreCase(n.getName())) {
			
			n.getParent().setRightChild(null);
			
			if(n.getParent().getLeftChild()!=null && !n.isWGD && !n.isWGT) {
				
				n.getParent().isVirtualNode=true; //for calculations of the likelihood

			}
			
			if(n.getParent().getLeftChild()==null && n.getParent().getRightChild()==null) {
				deleteNode(n.getParent());

			}
		}

	}

	public void deleteNodeFromTree(List<String> nodeNamesToDelete) {

		for (String s : nodeNamesToDelete) {
			
			this.deleteNode(this.findNodeWithName(s));

		}
	}

	// public void insertWGMsToTheTree_customMultFactor(List<List<String>> wgdList,
	// double[] retentionRates)
	//
	// {
	//
	// if(retentionRates.length != wgdList.size()){
	// System.out.println("WARNING: the length of multiplication factor's array and
	// number of WGDs do not match!");
	// }
	//
	// /* for each line of the list = each WGD:*/
	//
	// for (int i = 0; i < wgdList.size() ; i++) {
	//
	// String parentName = wgdList.get(i).get(1);
	// String childName = wgdList.get(i).get(2);
	//
	// Node wg = new Node();
	// wg.isWGM = true;
	// wg.setbLen(Double.parseDouble(wgdList.get(i).get(3)));
	// // becasue we call it with root.insertWGM.. , so the maxNodeSize of root is
	// applied, which is set in newickparser while building tree
	// wg.setmaxNodeSize(this.getmaxNodeSize());
	//
	// if (wgdList.get(i).get(0).equalsIgnoreCase("WGD")) {
	// wg.setName("WGD-" + parentName + "-" + childName);
	// wg.isWGD = true;
	// wg.multFactor = 1+(1*retentionRates[i]);
	// }
	//
	// if (wgdList.get(i).get(0).equalsIgnoreCase("WGT")) {
	// wg.setName("WGT-" + parentName + "-" + childName);
	// wg.isWGT = true;
	// wg.multFactor = 1+(2*retentionRates[i]);
	// }
	// addWGMRec(this, parentName, childName, wg);
	// }
	// }

	public void addWGMOnaSpecificBranch(Branch br, double blen) {

		String parentName = br.getParent().getName();
		String childName = br.getChild().getName();

		System.out.println("adding WGD to" + br.getParent().getName() + "  " + br.getChild().getName());
		Node wg = new Node();
		wg.setName("WGD " + parentName + "-" + childName);
		wg.isWGM = true;
		/* For now we only be.ugent.psb.setas.bdmodel.test for WGDs */
		wg.multFactor = 2;
		wg.setbLen(blen);
		wg.setmaxNodeSize(100);

		addDepthSubTree(br.getChild(), 1);
		wg.depth = (br.getParent().depth) + 1;

		double oldBlen = br.getChild().getbLen();
		br.getChild().setParent(wg);

		if (br.getParent().getLeftChild() != null) {

			if (br.getParent().getLeftChild().equals(br.getChild())) {

				wg.setLeftChild(br.getChild());
				wg.getLeftChild().setbLen(oldBlen - wg.getbLen());

				wg.setParent(br.getParent());
				br.getParent().setLeftChild(wg);
			}
		}

		else if (br.getParent().getRightChild() != null) {
			if (br.getParent().getRightChild().equals(br.getChild())) {

				wg.setRightChild(br.getChild());
				wg.getRightChild().setbLen(oldBlen - wg.getbLen());

				wg.setParent(br.getParent());
				br.getParent().setRightChild(wg);
			}
		}

		addWGMRec(this, parentName, childName, wg);

	}

	/* finds maximum element of a queue */
	public double findMaxQ(Queue<Double> q) {
		int size = q.size();
		double[] b = new double[size];

		for (int k = 0; k < size; k++) {
			b[k] = q.remove();
		}
		return FindMaxArray.findMaxValue(b);
	}

	public Queue<Double> allBlen() {
		Queue<Double> queue = new LinkedList<Double>();
		return allBlenRec(this, queue);
	}

	/** make a queue of all branch lengths in a tree */
	private Queue<Double> allBlenRec(Node n, Queue<Double> queue) {

		Node leftChild = n.getLeftChild();
		Node rightChild = n.getRightChild();

		if (leftChild != null) {
			queue.add(leftChild.getbLen());
			allBlenRec(leftChild, queue);
		}

		if (rightChild != null) {
			queue.add(rightChild.getbLen());
			allBlenRec(rightChild, queue);
		}

		return queue;
	}

	/** returns the maximum branch length of a tree rooted at Node n */
	public double maxBlen() {
		Queue<Double> qd = new LinkedList<Double>();
		qd = this.allBlen();
		double d = findMaxQ(qd);
		return d;
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

	public double calcMaxDisToRoot() {

		ArrayList<Node> arln = postOrder(this);

		double maxDist = arln.get(0).getDistToRoot();

		for (int i = 1; i < arln.size(); i++) {

			if (arln.get(i).getDistToRoot() > maxDist)

			{
				maxDist = arln.get(i).getDistToRoot();
			}
		}

		return maxDist;
	}

	/*
	 * To make all the leaves have the same distance to the root, for testing the
	 * program with comparing the results with CAFE
	 */
	public void adjustBranchLengths(ArrayList<Node> leaves) {

		double maxDisToRoot = this.calcMaxDisToRoot();

		for (int i = 0; i < leaves.size(); i++) {

			double oldBlen = leaves.get(i).getbLen();

			double distance = leaves.get(i).calcDistToRoot();

			double difference = maxDisToRoot - distance;

			double newBlen = oldBlen + difference;
			leaves.get(i).setbLen(newBlen);

		}
	}

	// public void setMaxSizes(Node r, int maxSize) {
	//
	// int sizeLkArray=0;
	// if (r != null) {
	// r.setmaxNodeSize(maxSize);
	//
	// if(r.isWGM==false){
	// sizeLkArray= maxSize;
	// }
	//
	// else{
	// sizeLkArray= (r.multFactor*maxSize);
	// }
	//
	// r.lk= new double[sizeLkArray+1];
	// Node left = r.getLeftChild();
	// Node right = r.getRightChild();
	// setMaxSizes(left,sizeLkArray);
	// setMaxSizes(right,sizeLkArray);
	//
	// }
	// }

	// same order as written in Newick format: left-right-root
	public ArrayList<Node> postOrder(Node root) {

		Stack<Node> stack = new Stack<Node>();
		Stack<Node> stackTwo = new Stack<Node>();

		ArrayList<Node> arln = new ArrayList<Node>();

		if (root != null) {
			stack.push(root);

			while (!stack.empty()) {
				root = stack.pop();
				stackTwo.push(root);

				if (root.getLeftChild() != null) {
					stack.push(root.getLeftChild());
				}
				if (root.getRightChild() != null) {
					stack.push(root.getRightChild());
				}
			}
		}

		while (!stackTwo.empty()) {
			arln.add(stackTwo.pop());
		}
		return arln;
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

	// public static void main(String args[]) {
	// NewickParser np = new NewickParser();
	// Node root = np
	// .buildTree("/home/setas/workspace/BirthDeathModel/src/Files/Trees/Pruned37To14ToEudicotTreeMGCF5");
	// ArrayList<Node> leaves = new ArrayList<Node>();
	// leaves = root.getLeaves();
	//
	// root.setNumberOfLeaves(leaves.size());
	//
	// int[] testObs = new int[root.numberOfLeaves];
	//
	// for (int i = 0; i < testObs.length; i++) {
	// testObs[i] = i + 1;
	// // gives the number of leaf because we start at 0
	// }
	// root.setLeafValues(testObs);
	//
	// WGMparser wgm = new WGMparser();
	// List<List<String>> wgdList = wgm
	// .readInputFile("/home/setas/workspace/BirthDeathModel/src/Files/WGMs/WGDs-avg-pruned37MGCF5ToEudicots-Clean");
	//
	// root.insertWGMsToTheTree(wgdList);
	//
	// ArrayList<Node> arln= root.postOrder(root);
	//
	//// for(int i=0; i<arln.size();i++){
	////
	//// System.out.println(arln.get(i).getName());
	//// }
	// double max = root.maxBlen();
	//
	//// To prepare the tree to be.ugent.psb.setas.bdmodel.test with cafe
	//// root.adjustBranchLengths(leaves);
	////
	//// ArrayList<Node> arp = root.postOrder(root);
	//// for(int i=0; i<arp.size();i++){
	////
	////// arp.get(i).setbLen(arp.get(i).getbLen()*1000000);
	////
	//// System.out.println(arp.get(i).getName()+" "+ arp.get(i).getbLen());
	//// }
	// ArrayList<Branch> arBranches= root.findAllBranches(root);
	//// System.out.println(arBranches.size());
	//
	//// for(int i=0; i<arBranches.size();i++){
	////
	//// System.out.println(arBranches.get(i).getParent()+ "
	// "+arBranches.get(i).getChild()+" "+ arBranches.get(i).getBranchLenght());ss
	//// }
	//
	// }
}