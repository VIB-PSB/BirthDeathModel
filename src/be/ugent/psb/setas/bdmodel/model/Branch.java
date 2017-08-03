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

/**
 * Defines all properties of a Branch Object, namely: Parent Node, Child Node,
 * Number of partitions on that branch Branch length that is left after making
 * equi-distance partitioning and inserting virtual nodes at those partitions A
 * boolean - isSelected - that is used to traverse all branches of the tree to
 * insert virtual nodes.
 */
public class Branch {


	private Node parent;
	private Node child;
	private double branchLength;
	private boolean isSelected;
	private int numberOfWGDs;

	private int numberOfPartitions;

	@Override
	public String toString() {
		return "Branch: " + this.getParent().getName() + "\t" + this.getChild().getName();
	}

	public Branch(Node parent, Node child) {
		this.parent = parent;
		this.child = child;
		this.branchLength = child.getbranchLength();
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
		this.branchLength = this.child.getbranchLength();
	}

	public boolean getIsSelected() {
		return isSelected;
	}

	public void setIsSelected(boolean b) {
		this.isSelected = b;
	}

}
