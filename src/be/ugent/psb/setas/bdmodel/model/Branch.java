package be.ugent.psb.setas.bdmodel.model;

public class Branch {

	private Node parent;
	private Node child;
	private double branchLength;
	private boolean isSelected;
	private double mostLikelyPlaceWGD;
	private double[] mostLikelyPlacesOf2WGDs;
	private int numberOfWGDs ;
	
	private int numberOfPartitions;
	
	
	@Override
	public String toString() {
		return "Branch: " + this.getParent().getName()+"\t"+this.getChild().getName();
	}
	
	
	public Branch(Node parent, Node child) {
		this.parent = parent;
		this.child = child;
		this.branchLength = child.getbLen();
	}
	
	
	public int getNumberOfPartitions() {
		return numberOfPartitions;
	}

	public void setNumberOfPartitions(int numberOfPartitions) {
		this.numberOfPartitions = numberOfPartitions;
	}

	private double remainingBlen;
	
	
	
	public double getRemainingBlen() {
		return remainingBlen;
	}

	public void setRemainingBlen(double remainingBlen) {
		this.remainingBlen = remainingBlen;
	}

	public int getNumberOfWGDs() {
		return numberOfWGDs;
	}

	public void setNumberOfWGDs(int numberOfWGDs) {
		this.numberOfWGDs = numberOfWGDs;
	}

	public double[] getMostLikelyPlacesOf2WGDs() {
		return mostLikelyPlacesOf2WGDs;
	}

	public void setMostLikelyPlacesOf2WGDs(double[] mostLikelyPlacesOf2WGDs) {
		this.mostLikelyPlacesOf2WGDs = mostLikelyPlacesOf2WGDs;
	}

	public double getMostLikelyPlaceWGD() {
		return mostLikelyPlaceWGD;
	}

	public void setMostLikelyPlaceWGD(double mostLikelyPlaceWGD) {
		this.mostLikelyPlaceWGD = mostLikelyPlaceWGD;
	}


	public Node getParent() {
		return parent;
	}

	public void setParent(Node parent) {
		this.parent = parent;
	}

	public Node getChild() {
		return child;
	}

	public void setChild(Node c) {
		this.child = c;
	}
	
	public double getBranchLenght() {
		return branchLength;
	}

	public void setBranchLength() {
		this.branchLength = this.child.getbLen();
	}
	
	public boolean getIsSelected() {
		return isSelected;
	}

	public void setIsSelected(boolean b) {
		this.isSelected = b;
	}
	
}
