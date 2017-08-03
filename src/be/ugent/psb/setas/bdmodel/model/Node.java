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
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Stack;

public class Node {

	public Node() {
	}

	private Node leftChild;
	private Node rightChild;
	private Node parentNode;

	/**
	 * If a node has no children, it is called a leaf and the boolean isleaf is
	 * True for it.
	 */
	public boolean isRoot;
	public boolean isLeaf;
	public boolean isWGD;
	public boolean isWGT;
	public boolean isWGM;
	public boolean isVirtualNode;

	/**
	 * multiplication factor is equal to 2 for WGDs and equal to 3 for WGMs. It is not defined or used for other nodes.
	 */
	public int multiplicationFactor;

	private String nodeName;

	/** The distance between a node and its parentNode */
	private double branchLength;
	private double distanceToRoot;

	private int geneCountAtNode;

	public int getGeneCountAtNode() {
		return this.geneCountAtNode;
	}

	public void setValue(int value) {
		this.geneCountAtNode = value;
	}

	public int depthOfNode;

	private int maximumGeneCountAtNode;
	private double[] likelihoodArrayOfNode = new double[maximumGeneCountAtNode + 1];

	private ArrayList<Node> arrayListOfLeaves = new ArrayList<Node>();

	/**
	 * The array reference likelihood saves likelihoods of the geneCountAtNode of a node
	 * being 1 to maximumGeneCount at nodes (set to 100 in the main method), 
	 * given an array of observations from the input file (= reference observation).
	 */
	public double[] getLikelihoodArray() {

		return this.likelihoodArrayOfNode;
	}

	public double getLikelihoodValue(int geneCountAtTheNode) {

		return this.likelihoodArrayOfNode[geneCountAtTheNode];
	}

	public void setLikelihoodValue(int geneCount, double likelihood) {
		if (geneCount > this.likelihoodArrayOfNode.length) {
			geneCount = this.maximumGeneCountAtNode;
		}
		likelihoodArrayOfNode[geneCount] = likelihood;
	}

	public int getMaxGeneCountAtNode() {
		return maximumGeneCountAtNode;
	}

	public void setMaxNodeGeneCountAtNode(int maximumGeneCount) {
		this.maximumGeneCountAtNode = maximumGeneCount;
		double[] newLogklikelihood = new double[this.maximumGeneCountAtNode + 1];
		for (int i = 0; i < likelihoodArrayOfNode.length; i++) {
			newLogklikelihood[i] = likelihoodArrayOfNode[i];
		}
		this.likelihoodArrayOfNode = newLogklikelihood;
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
		return parentNode;
	}

	public void setParent(Node parent) {
		this.parentNode = parent;
	}

	public String getName() {
		return nodeName;
	}

	public void setName(String nameOfNode) {
		this.nodeName = nameOfNode;
	}


	/**
	 * @return the branch length between a node and its parentNode
	 */
	public double getbranchLength() {
		return branchLength;
	}

	public void setbranchLength(double branchlength) {
		this.branchLength = branchlength;
	}

	public double getDistanceToRoot() {
		return distanceToRoot;
	}

	public void setDistanceToRoot(double distanceToRoot) {

		this.distanceToRoot = distanceToRoot;
	}

	@Override
	public String toString() {
		return "Node: " + this.getName();
	}

	public int getNumberOfLeaves() {
		return arrayListOfLeaves.size();
	}

	/**
	 * @return list of leaf nodes in post order:  left-right-(root) 
	 */
	public ArrayList<Node> getLeaves() {
		Stack<Node> tempStack = new Stack<Node>();
		Stack<Node> stackOfLeaves = new Stack<Node>();

		if (this != null) {
			tempStack.push(this);

			while (!tempStack.empty()) {
				Node root = tempStack.pop();

				if (root.getLeftChild() == null && root.getRightChild() == null) {
					stackOfLeaves.push(root);
				}

				if (root.getLeftChild() != null) {
					tempStack.push(root.getLeftChild());
				}
				if (root.getRightChild() != null) {
					tempStack.push(root.getRightChild());
				}
			}
		}

		arrayListOfLeaves = new ArrayList<Node>();
		while (!stackOfLeaves.empty()) {
			this.arrayListOfLeaves.add(stackOfLeaves.pop());
		}
		return this.arrayListOfLeaves;
	}

	/**
	 * CAUTION: The species in gene family counts' file should be in the same order as in newick format
	 * file, i.e. post order
	 */
	public void setLeafValues(int[] geneCountArray) {
		if (arrayListOfLeaves.isEmpty()) {
			this.getLeaves();
		}

		if (geneCountArray.length == this.getNumberOfLeaves()) {
			for (int i = 0; i < geneCountArray.length; i++) {
				this.arrayListOfLeaves.get(i).setValue(geneCountArray[i]);

			}
		} else {
			System.err.println("Error: GF-Counts' length does not matche the number of leaves!");
		}
	}

	public void addDepthSubTree(Node node, int addToTheDepth) {

		if (node != null) {
			node.depthOfNode += addToTheDepth;
			addDepthSubTree(node.getLeftChild(), addToTheDepth);
			addDepthSubTree(node.getRightChild(), addToTheDepth);
		}
	}

	/**
	 * 
	 * @return All nodes in the tree with this node as root, traversed in
	 *         postorder
	 */
	public Queue<Node> postOrder() {
		Queue<Node> queueOfNodes = new LinkedList<Node>();
		return recurviseTraverseTreeInPostOrder(this, queueOfNodes);
	}

	private Queue<Node> recurviseTraverseTreeInPostOrder(Node root, Queue<Node> queueOfNodesInPostOrder) {
		if (root != null) {
			recurviseTraverseTreeInPostOrder(root.getLeftChild(), queueOfNodesInPostOrder);
			recurviseTraverseTreeInPostOrder(root.getRightChild(), queueOfNodesInPostOrder);
			queueOfNodesInPostOrder.add(root);
		}
		return queueOfNodesInPostOrder;
	}

	public Node findNodeWithName(String nameOfNodeToSearch) {
		Queue<Node> allNodes = this.postOrder();
		for (Node node : allNodes) {
			if (node.getName().equals(nameOfNodeToSearch)) {
				return node;
			}
		}
		return null;
	}

	public void addWGM(String parentName, String childName, Node wgmNode) {
		addWGMRecursively(this, parentName, childName, wgmNode);
	}

	/**
	 * If there are subsequent WGMs, in order for this function to work
	 * correctly, they should be mentioned in the data file in order from oldest
	 * to the youngest.
	 */
	private void addWGMRecursively(Node node, String parentsName, String childsName, Node wgmNode) {

		if (node != null) {

			if (!node.isLeaf) {

				if (node.getName().equals(parentsName)) {

					if (node.getLeftChild() != null && node.getLeftChild().getName().equals(childsName)) {

						addDepthSubTree(node.getLeftChild(), 1);
						wgmNode.depthOfNode = (node.depthOfNode) + 1;

						double oldBlen = node.getLeftChild().getbranchLength();

						node.getLeftChild().setParent(wgmNode);

						wgmNode.setLeftChild(node.getLeftChild());
						wgmNode.getLeftChild().setbranchLength(oldBlen - wgmNode.getbranchLength());

						wgmNode.setParent(node);
						node.setLeftChild(wgmNode);
					}

					if (node.getRightChild() != null && node.getRightChild().getName().equals(childsName)) {

						addDepthSubTree(node.getRightChild(), 1);
						wgmNode.depthOfNode = (node.depthOfNode) + 1;

						double oldBlen = node.getRightChild().getbranchLength();

						node.getRightChild().setParent(wgmNode);
						wgmNode.setRightChild(node.getRightChild());
						wgmNode.getRightChild().setbranchLength(oldBlen - wgmNode.getbranchLength());

						wgmNode.setParent(node);
						node.setRightChild(wgmNode);
					}
				}

				else {
					Node left = node.getLeftChild();
					Node right = node.getRightChild();
					if (left != null) {
						addWGMRecursively(left, parentsName, childsName, wgmNode);
					}
					if (right != null) {
						addWGMRecursively(right, parentsName, childsName, wgmNode);
					}
				}
			}
		}

	}

	/**
	 * To read wgds from a list and build the corresponding nodes and use
	 * addWGMRec method to add them recursively to the tree
	 */
	public void insertWGMsToTheTree(List<List<String>> wgmList)

	{
		for (int i = 0; i < wgmList.size(); i++) {

			String parentName = wgmList.get(i).get(1);
			String childName = wgmList.get(i).get(2);

			Node wgm = new Node();
			wgm.isWGM = true;
			wgm.setbranchLength(Double.parseDouble(wgmList.get(i).get(3)));
			wgm.setMaxNodeGeneCountAtNode(this.getMaxGeneCountAtNode());

			if (wgmList.get(i).get(0).equalsIgnoreCase("WGD")) {
				wgm.setName("WGD-" + parentName + "-" + childName);
				wgm.isWGD = true;
				wgm.multiplicationFactor = 2;
			}

			if (wgmList.get(i).get(0).equalsIgnoreCase("WGT")) {
				wgm.setName("WGT-" + parentName + "-" + childName);
				wgm.isWGT = true;
				wgm.multiplicationFactor = 3;
			}
			addWGMRecursively(this, parentName, childName, wgm);
		}
	}
	
	/**
	 * Returns the length of path from a lead node to the root
	 * @return
	 */

	public double calculateDistanceToRoot() {

		double distance = this.getbranchLength();

		Node parentNode = this.getParent();

		while (!parentNode.isRoot) {

			distance += parentNode.getbranchLength();

			parentNode = parentNode.getParent();
		}
		return distance;
	}

	/** same order as written in Newick format: left-right-root **/
	public ArrayList<Node> postOrder(Node root) {

		Stack<Node> stackOneTemp = new Stack<Node>();
		Stack<Node> stackTwoTemp = new Stack<Node>();

		ArrayList<Node> arrayListOfNodesInPostOrder = new ArrayList<Node>();

		if (root != null) {
			stackOneTemp.push(root);

			while (!stackOneTemp.empty()) {
				root = stackOneTemp.pop();
				stackTwoTemp.push(root);

				if (root.getLeftChild() != null) {
					stackOneTemp.push(root.getLeftChild());
				}
				if (root.getRightChild() != null) {
					stackOneTemp.push(root.getRightChild());
				}
			}
		}

		while (!stackTwoTemp.empty()) {
			arrayListOfNodesInPostOrder.add(stackTwoTemp.pop());
		}
		return arrayListOfNodesInPostOrder;
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
