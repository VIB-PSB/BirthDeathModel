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
import java.util.Stack;

public class PartitionBranches {

	public PartitionBranches(Node root, double desiredmaximumLenghtOfBranches) {
		this.root = root;
		this.maximumLenghtOfBranches = desiredmaximumLenghtOfBranches;
	}

	private Node root;
	private double maximumLenghtOfBranches;

	public double getThreshold() {
		return maximumLenghtOfBranches;
	}

	public void setThreshold(double threshold) {
		this.maximumLenghtOfBranches = threshold;
	}

	public Node getRoot() {
		return root;
	}

	public void setRoot(Node root) {
		this.root = root;
	}

	public void setParametersForBranchLongerThanThreshold(Branch branch) {

		double branchLength = branch.getBranchLenght();

		if (branchLength > maximumLenghtOfBranches) {

			int numberOfPartitions = (int) (Math.floor(branchLength / maximumLenghtOfBranches));
			double remainingBranchlength = branchLength - (numberOfPartitions) * maximumLenghtOfBranches;

			branch.setNumberOfPartitions(numberOfPartitions);
			branch.setRemainingBlen(remainingBranchlength);
		}

	}

	public Node addVirtualNodeOnaSpecificBranch(Branch branch) {

		setParametersForBranchLongerThanThreshold(branch);

		Node child = branch.getChild();
		Node parent = branch.getParent();

		String parentName = parent.getName();
		String childName = child.getName();

		double oldBlen = child.getbranchLength();

		Node virtualNode = new Node();
		virtualNode.isVirtualNode = true;
		virtualNode.setName("VN " + parentName + " to " + childName);
		virtualNode.setbranchLength(maximumLenghtOfBranches);
		virtualNode.depthOfNode = (parent.depthOfNode) + 1;
		virtualNode.setMaxNodeGeneCountAtNode(parent.getMaxGeneCountAtNode());

		child.addDepthSubTree(child, 1);
		child.setParent(virtualNode);

		if (parent.getLeftChild() != null
				&& parent.getLeftChild().getName().equals(child.getName())) {

			virtualNode.setLeftChild(child);
			virtualNode.getLeftChild().setbranchLength(oldBlen - virtualNode.getbranchLength());

			virtualNode.setParent(parent);
			parent.setLeftChild(virtualNode);
		}

		else if (parent.getRightChild() != null
				&& parent.getRightChild().getName().equals(child.getName())) {

			virtualNode.setRightChild(child);
			virtualNode.getRightChild().setbranchLength(oldBlen - virtualNode.getbranchLength());

			virtualNode.setParent(parent);
			parent.setRightChild(virtualNode);
		}

		return virtualNode;

	}

	public void addAllVirtualNodesOnABranch(Branch branch) {

		Stack<Branch> stackOfBranches = new Stack<Branch>();
		addAllVirtualNodesOnABranchRecursively(branch, stackOfBranches);
	}

	private void addAllVirtualNodesOnABranchRecursively(Branch newBranch, Stack<Branch> stackOfBranches) {

		stackOfBranches.push(newBranch);

		while (!stackOfBranches.empty()) {

			Branch testBranch = stackOfBranches.pop();

			if (testBranch.getNumberOfPartitions() > 0) {

				Node virtualNode = addVirtualNodeOnaSpecificBranch(testBranch);

				if (virtualNode.getLeftChild() != null) {
					Branch leftBranch = new Branch(virtualNode, virtualNode.getLeftChild());
					
					setParametersForBranchLongerThanThreshold(leftBranch);					
					stackOfBranches.push(leftBranch);

				}

				else if (virtualNode.getRightChild() != null) {

					Branch rightBranch = new Branch(virtualNode, virtualNode.getRightChild());
					setParametersForBranchLongerThanThreshold(rightBranch);
					stackOfBranches.push(rightBranch);
				}
			}
		}

	}

	public void insertAllVirtualNodesOnAllBranches() {

		ArrayList<Branch> arrayListOfBranches = root.findAllBranches(root);

		for (Branch branch : arrayListOfBranches) {

			setParametersForBranchLongerThanThreshold(branch);
			addAllVirtualNodesOnABranch(branch);
		}
	}

}
