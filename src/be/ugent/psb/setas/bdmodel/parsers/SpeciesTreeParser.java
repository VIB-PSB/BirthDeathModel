package be.ugent.psb.setas.bdmodel.parsers;

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
import java.util.List;
import java.util.Queue;

import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.bdmodel.model.PartitionBranches;

public class SpeciesTreeParser {

	/**
	 * Given the newick format tree, WGMs and partition size of branches, 
	 * This class builds the tree, inserts the WGMs 
	 * and partition the branches of the tree 
	 * and inserts Virtual Nodes at those partitions.
	 * @param treeFile
	 * @param wgmFile
	 * @param partitionSize
	 * @param maximumGeneCount
	 * @return
	 */
	public static Node buildAndPartitionTree(String treeFile, String wgmFile, double partitionSize,
			int maximumGeneCount) {

		NewickFormatTreeParser newickParser = new NewickFormatTreeParser();
		Node root = newickParser.setParametersOfTree(treeFile, maximumGeneCount);
		root.getLeaves();

		WGMparser wgm = new WGMparser();
		List<List<String>> wgmList = wgm.readWGMfile(wgmFile);
		root.insertWGMsToTheTree(wgmList);

		if (partitionSize != 0) {
			PartitionBranches partitionBranches = new PartitionBranches(root, partitionSize);
			partitionBranches.insertAllVirtualNodesOnAllBranches();
		}

		return root;
	}
	
	/**
	 * Sets the values of gene counts for the spcified gene family
	 * @param root
	 * @param allGeneCounts
	 * @param numberOfGeneFamily
	 */
	public static void setLeavesValues(Node root, List<List<Integer>> allGeneCounts, int numberOfGeneFamily) {

		List<Integer> listOfGeneCounts = allGeneCounts.get(numberOfGeneFamily);

		int[] geneCountProfile = new int[listOfGeneCounts.size()];

		for (int speciesNumber = 0; speciesNumber < listOfGeneCounts.size(); speciesNumber++) {
			
			geneCountProfile[speciesNumber] = listOfGeneCounts.get(speciesNumber);
		}
		root.setLeafValues(geneCountProfile);
	}

	public static ArrayList<Node> setMaximumGeneCount(Node root, int maximumGeneCount) {

		Queue<Node> queueOfNodesInPostOrder = root.postOrder();

		for (Node node : queueOfNodesInPostOrder) {
			node.setMaxNodeGeneCountAtNode(maximumGeneCount);
		}

		ArrayList<Node> arrayListOfNodes = new ArrayList<Node>(queueOfNodesInPostOrder);
		return arrayListOfNodes;

	}
}
