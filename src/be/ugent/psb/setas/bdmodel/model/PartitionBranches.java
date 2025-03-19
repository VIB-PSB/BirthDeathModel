package be.ugent.psb.setas.bdmodel.model;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Stack;

import be.ugent.psb.setas.bdmodel.parsers.NewickParser;

public class PartitionBranches {

	public PartitionBranches(Node root, double threshold) {
		this.root = root;
		this.threshold = threshold;
	}

	private Node root;
	private double threshold; // max t in the tree we wish to use, e.g. 0.1

	public double getThreshold() {
		return threshold;
	}

	public void setThreshold(double threshold) {
		this.threshold = threshold;
	}

	public Node getRoot() {
		return root;
	}

	public void setRoot(Node root) {
		this.root = root;
	}

	public void setParamsForBrLongerThanThrshold(Branch br) {

		double blen = br.getBranchLenght();

		if (blen > threshold) {

			int numOfPartitions = (int) (Math.floor(blen / threshold));
			double remainingBlen = blen - (numOfPartitions) * threshold;

			br.setNumberOfPartitions(numOfPartitions);
			br.setRemainingBlen(remainingBlen);
		}

	}

	public Node addVirtualNodeOnaSpecificBranch(Branch br) {

		setParamsForBrLongerThanThrshold(br);

		Node child = br.getChild();
		Node parent = br.getParent();

		String parentName = parent.getName();
		String childName = child.getName();

		double oldBlen = child.getbLen();

		Node vn = new Node();
		vn.isVirtualNode = true;
		vn.setName("VN " + parentName + " to " + childName);
		vn.setbLen(threshold);
		vn.depth = (parent.depth) + 1;
		vn.setmaxNodeSize(parent.getmaxNodeSize());

		child.addDepthSubTree(child, 1);
		child.setParent(vn);

		if (parent.getLeftChild() != null
				&& parent.getLeftChild().getName().equals(child.getName())) {

			vn.setLeftChild(child);
			vn.getLeftChild().setbLen(oldBlen - vn.getbLen());

			vn.setParent(parent);
			parent.setLeftChild(vn);
		}

		else if (parent.getRightChild() != null
				&& parent.getRightChild().getName().equals(child.getName())) {

			vn.setRightChild(child);
			vn.getRightChild().setbLen(oldBlen - vn.getbLen());

			vn.setParent(parent);
			parent.setRightChild(vn);
		}

		return vn;

	}

	public void addAllVirtualNodesOnABranch(Branch br) {

		Stack<Branch> stackOfBranches = new Stack<Branch>();
		addAllVirtualNodesOnABranchRec(br, stackOfBranches);
	}

	private void addAllVirtualNodesOnABranchRec(Branch br, Stack<Branch> stBr) {

		stBr.push(br);

		while (!stBr.empty()) {

			Branch test = stBr.pop();

			if (test.getNumberOfPartitions() > 0) {

				Node vn = addVirtualNodeOnaSpecificBranch(test);

				if (vn.getLeftChild() != null) {
					Branch brLeft = new Branch(vn, vn.getLeftChild());
					
					setParamsForBrLongerThanThrshold(brLeft);					
					stBr.push(brLeft);

				}

				else if (vn.getRightChild() != null) {

					Branch brRight = new Branch(vn, vn.getRightChild());
					setParamsForBrLongerThanThrshold(brRight);
					stBr.push(brRight);
				}
			}
		}

	}

	public void insertAllVNSonAllBranches() {

		ArrayList<Branch> arb = root.findAllBranches(root);

		for (Branch br : arb) {

			setParamsForBrLongerThanThrshold(br);
			addAllVirtualNodesOnABranch(br);
		}
	}

//	 public static void main (String [] args){
//	
//	 int maxNodeSize = 5;
//	
//	 NewickParser np = new NewickParser();
//	 Node root = np.buildTree(args[0],maxNodeSize);
//	
//	 ArrayList<Node> leaves = root.getLeaves();
//	 root.setLeaves(leaves);
//	 root.setNumberOfLeaves(leaves.size());
//	
//	 PartitionBranches pb = new PartitionBranches(root,0.1);
//	
//	 pb.insertAllVNSonAllBranches();
//	
//	 Queue<Node> queue = root.postOrder();
//	 ArrayList<Node> aryln = new ArrayList<Node>(queue);
//	
//	 for(Node n : aryln){
//	
//	 System.out.println(n.getName());
//	 }
//	
//	
//	
//	 }
}
