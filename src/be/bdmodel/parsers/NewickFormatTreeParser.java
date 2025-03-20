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

import java.io.FileReader;
import java.util.ArrayList;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Scanner;
import java.util.Stack;

import be.ugent.psb.setas.bdmodel.model.Node;

public class NewickFormatTreeParser {

	public String[] queueOfStringsToArrayOfString(Queue<String> queueOfStrings) {
		int size = queueOfStrings.size();
		String[] arrayOfStrings = new String[size];

		for (int k = 0; k < size; k++) {
			arrayOfStrings[k] = queueOfStrings.remove();
		}
		return arrayOfStrings;
	}

	/**
	 * reads the newick input tree as a string and removes all ( and ) and ; and / and :
	 * @param newickFormatFile
	 * @return
	 */
	public String[] readNewickFormatFile(String newickFormatFile) {

		FileReader fileReader = null;
		try {
			fileReader = new FileReader(newickFormatFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner scanner = new Scanner(fileReader);
		String line = new String();

		Queue<String> queueOfLineParts = new LinkedList<String>();

		while (scanner.hasNextLine()) {
			line = scanner.nextLine();
			String[] parts = line.split("(\\,)|(:)|(\\()|(\\))|(;)");

			for (int i = 0; i < parts.length; i++) {
				queueOfLineParts.add(parts[i]);
			}
		}
		scanner.close();

		String[] nodesAndBranchLengthsArray = queueOfStringsToArrayOfString(queueOfLineParts);

		return nodesAndBranchLengthsArray;
	}

	public ArrayList<Node> buildNodes(String[] nodesAndBranchLengthsArray) {

		int numberOfNodesAndBranchLengths = nodesAndBranchLengthsArray.length;

		Queue<String> queueOfNamesAndBranchLenghts = new LinkedList<String>();
		ArrayList<Node> arrayListOfNodes = new ArrayList<Node>();

		for (int k = 0; k < numberOfNodesAndBranchLengths; k++) {
			if (!nodesAndBranchLengthsArray[k].isEmpty()) {
				queueOfNamesAndBranchLenghts.add(nodesAndBranchLengthsArray[k]);
			}
		}

		String[] nodeInformation = queueOfStringsToArrayOfString(queueOfNamesAndBranchLenghts);

		int lenxOfNodeInformation = nodeInformation.length;
		for (int i = 0; i < lenxOfNodeInformation - 1; i = i + 2) {
			Node node = new Node();
			node.setName(nodeInformation[i]);
			node.setbranchLength(Double.parseDouble(nodeInformation[i + 1]));
			arrayListOfNodes.add(node);
					}
		Node root = new Node();
		root.setName(nodeInformation[lenxOfNodeInformation - 1]);
		arrayListOfNodes.add(root);
		return arrayListOfNodes;
	}

	public static boolean isNumeric(String string) {
		return string.matches("-?\\d+(\\.\\d+)?"); // match a number with optional
												// '-' and decimal.
	}

	/** WITH parenthesis, without numbers: (without , and : and ; )
	 * 
	 * @param newickFormatTree
	 * @return
	 */
	public String[] readNewickTree_KeepParanthesis(String newickFormatTree) {

		FileReader fileReader = null;
		try {
			fileReader = new FileReader(newickFormatTree);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner scanner = new Scanner(fileReader);
		String newickString = new String();

		Queue<String> queueOfNodeNameIncludingParanthsis = new LinkedList<String>();

		while (scanner.hasNextLine()) {
			
			newickString = scanner.nextLine();
			String[] parts = newickString.split("(\\,)|(:)|(;)");

			for (int i = 0; i < parts.length; i++) {
				if (!isNumeric(parts[i])) {
					queueOfNodeNameIncludingParanthsis.add(parts[i]);
				}
			}
		}
		scanner.close();

		String[] nodeNamesIncludingParanthesis = queueOfStringsToArrayOfString(queueOfNodeNameIncludingParanthsis);

		return nodeNamesIncludingParanthesis;
	}

	public int[] calculateDepth(String[] nodeNamesIncludingParanthesis) {

		int lenx = nodeNamesIncludingParanthesis.length;
		int[] depthOfNodes = new int[lenx];

		int nodeDepth = 0;

		for (int i = 0; i < lenx; i++) {
			
			String oneNodeNameIncludingParantheses = nodeNamesIncludingParanthesis[i];
			
			int numberOfCharactersInNodeNameIncludingParantheses = oneNodeNameIncludingParantheses.length();
			char[] charactersInNodeName = new char[numberOfCharactersInNodeNameIncludingParantheses];

			for (int j = 0; j < numberOfCharactersInNodeNameIncludingParantheses; j++) {

				charactersInNodeName[j] = oneNodeNameIncludingParantheses.charAt(j);
				if (charactersInNodeName[j] == '(') {
					nodeDepth += 1;
				}
				if (charactersInNodeName[j] == ')') {
					nodeDepth -= 1;
				}
			}
			depthOfNodes[i] = nodeDepth;
		}
		 for(int i=0;i<depthOfNodes.length;i++){
		 }
		return depthOfNodes;
	}

	public void setDepth(ArrayList<Node> arrayListOfNodes, String[] nodeNamesIncludingParanthesis) {

		int[] depths = calculateDepth(nodeNamesIncludingParanthesis);
		int numberOfNodes = arrayListOfNodes.size();
		
		for (int nodeN = 0; nodeN < numberOfNodes; nodeN++) {
			
			arrayListOfNodes.get(nodeN).depthOfNode = depths[nodeN];			
		}
	}

	
	/**
	 * Following the structure of newick format: first node: left most node and so on
	 * @param arratListOfNodes
	 * @return
	 */
	public Node buildTree(ArrayList<Node> arratListOfNodes) {
		Node root = new Node();
		Stack<Node> stackOfNodesTemp = new Stack<Node>();

		stackOfNodesTemp.push(arratListOfNodes.get(0));

		for (int i = 1; i < arratListOfNodes.size(); i++) {

			Node currentNode = arratListOfNodes.get(i);

			if (currentNode.depthOfNode >= stackOfNodesTemp.peek().depthOfNode) {
				stackOfNodesTemp.push(currentNode);
			}

			else {
				stackOfNodesTemp.peek().setParent(currentNode);
				currentNode.setRightChild(stackOfNodesTemp.pop());
				
				stackOfNodesTemp.peek().setParent(currentNode);
				currentNode.setLeftChild(stackOfNodesTemp.pop());
				stackOfNodesTemp.push(currentNode);
			}

		}

		root = stackOfNodesTemp.pop();
		root.isRoot =true;

		return root;
	}
	
	
	
	public void setDistanceToRoot(ArrayList<Node> arrayListOfNodes) {

		int size = arrayListOfNodes.size();
		
		for (int numberOfNode = 0; numberOfNode < size; numberOfNode++) {

			if(!arrayListOfNodes.get(numberOfNode).isRoot){
			
			double distanceToRoot= arrayListOfNodes.get(numberOfNode).calculateDistanceToRoot();
			
			arrayListOfNodes.get(numberOfNode).setDistanceToRoot(distanceToRoot);
			}
			
			else{
				arrayListOfNodes.get(numberOfNode).setDistanceToRoot(0);
			}
		}
	}

	/**
	 * Set Boolean isLeaf for all nodes
	 * @param node
	 */
	public void setIsleafBoolean(Node node) {

		if (node != null) {
			Node left = node.getLeftChild();
			Node right = node.getRightChild();

			if (right == null && left == null) {
				node.isLeaf = true;
			}

			setIsleafBoolean(left);
			setIsleafBoolean(right);
		}
	}
	
	public void setMaximumGeneCount(ArrayList<Node> arrayListOfNodes, int defaultMaximumGeneCount){
		
	for (int nodeNumber = 0; nodeNumber < arrayListOfNodes.size(); nodeNumber++) {
		arrayListOfNodes.get(nodeNumber).setMaxNodeGeneCountAtNode(defaultMaximumGeneCount);
	}
	}
	public Node setParametersOfTree(String newickFormatFile , int defaultMaximumGeneCount) {

		String[] nodesAndBranchLengths = readNewickFormatFile(newickFormatFile);
		String[] nodesAndBranchLengthsWithParanthesis = readNewickTree_KeepParanthesis(newickFormatFile);

		ArrayList<Node> arrayListOfNodes = buildNodes(nodesAndBranchLengths);
		setDepth(arrayListOfNodes, nodesAndBranchLengthsWithParanthesis);
	
		Node root = buildTree(arrayListOfNodes);
		
		setDistanceToRoot(arrayListOfNodes);			
		setMaximumGeneCount(arrayListOfNodes,defaultMaximumGeneCount);		
		setIsleafBoolean(root);
		
		return root;
	}


}
