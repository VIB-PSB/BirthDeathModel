package be.ugent.psb.setas.bdmodel.parsers;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import be.ugent.psb.setas.bdmodel.model.FindMaxArray;
import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.bdmodel.model.PartitionBranches;

public class SpeciesTreeParser {

	public static Node buildAndPartitionTree_NoWGD(String treeFile, double partitionSize, int defaultmaxNodeSize) {

		NewickParser np = new NewickParser();
		Node root = np.buildTree(treeFile, defaultmaxNodeSize);
		root.getLeaves();

		if (partitionSize != 0) {
			PartitionBranches pb = new PartitionBranches(root, partitionSize);
			pb.insertAllVNSonAllBranches();
		}
		return root;
	}

	public static Node buildInsertWGDsandPartitionTree(String treeFile, String wgdFile, double partitionSize,
			int defaultmaxNodeSize) {

		NewickParser np = new NewickParser();
		Node root = np.buildTree(treeFile, defaultmaxNodeSize);
		root.getLeaves();

		WGMparser wgm = new WGMparser();
		List<List<String>> wgdList = wgm.readInputFile(wgdFile);
		root.insertWGMsToTheTree(wgdList);

		if (partitionSize != 0) {
			PartitionBranches pb = new PartitionBranches(root, partitionSize);
			pb.insertAllVNSonAllBranches();
		}

		return root;
	}

//	public static Node  buildAndPartitionTree_customMultFactor(String treeFile, String wgdFile, String WGMretentionRatesFile, double partitionSize, int defaultmaxNodeSize) {
//
//		NewickParser np = new NewickParser();
//		Node root = np.buildTree(treeFile , defaultmaxNodeSize);
//		root.getLeaves();
//
//		WGMparser wgm = new WGMparser();
//		List<List<String>> wgdList = wgm.readInputFile(wgdFile);
//		ArrayList<String> retentionRates = wgm.readRetentionRates_String(WGMretentionRatesFile);
//		
//		double[] retRates = new double [retentionRates.size()];
//		
//		for(int i=0; i<retentionRates.size(); i++){
//			
//			retRates [i] = Double.parseDouble(retentionRates.get(i));
//		}
//		root.insertWGMsToTheTree_customMultFactor(wgdList, retRates); 
//		
//		if(partitionSize != 0){
//		PartitionBranches pb = new PartitionBranches(root, partitionSize);
//		pb.insertAllVNSonAllBranches();
//		}
//
//		return root;
//	}
//	
	public static Node buildAndPartitionTree_ReverseEng(String treeFile, String wgdFile, double partitionSize,
			int defaultmaxNodeSize, int ignoreLine) {

		NewickParser np = new NewickParser();
		Node root = np.buildTree(treeFile, defaultmaxNodeSize);
		root.getLeaves();

		WGMparser wgm = new WGMparser();
		List<List<String>> wgdList = wgm.readInputFile(wgdFile);
		root.insertWGMsToTheTree_ReverseEng(wgdList, ignoreLine);

		if (partitionSize != 0) {
			PartitionBranches pb = new PartitionBranches(root, partitionSize);
			pb.insertAllVNSonAllBranches();
		}
		return root;
	}

	public static Node buildAndPartitionTree_RVS(String treeFile, String wgdFileOrg, String rmWGD8, String rmWGD9,
			String rmWGD12, String rmWGD16, double partitionSize, int defaultmaxNodeSize, int numberOfWGDinList) {

		String wgdFileToPass = wgdFileOrg;

		if (numberOfWGDinList == 8) {wgdFileToPass = rmWGD8;}

		if (numberOfWGDinList == 9) {wgdFileToPass = rmWGD9;}

		if (numberOfWGDinList == 12) {wgdFileToPass = rmWGD12;}

		if (numberOfWGDinList == 16) {wgdFileToPass = rmWGD16;}

		return buildAndPartitionTree_ReverseEng(treeFile, wgdFileToPass, partitionSize, defaultmaxNodeSize,
				numberOfWGDinList);
	}

	public static void setLeavesValues(Node root, List<List<Integer>> nonRepetetiveCounts, int gf) {

		List<Integer> li = nonRepetetiveCounts.get(gf);
		// to exclude the gene counts of the last 9 non eudicot species of the 37 spe
//		int[] originalObservation = new int[li.size()-9];
		int[] originalObservation = new int[li.size()];

//		for (int m = 0; m < li.size()-9; m++) {
		for (int m = 0; m < li.size(); m++) {
			originalObservation[m] = li.get(m);
		}
		root.setLeafValues(originalObservation);
	}

	public static void setLeavesValues_one(Node root, List<Integer> li) {

		int[] originalObservation = new int[li.size()];

		for (int m = 0; m < li.size(); m++) {originalObservation[m] = li.get(m);}
		root.setLeafValues(originalObservation);
	}

	public static ArrayList<Node> setMaxNodeSizeAccToGF(Node root, List<List<Integer>> nonRepetetiveCounts, int gf) {

		Queue<Node> queue = root.postOrder();

		List<Integer> li = nonRepetetiveCounts.get(gf);
		int[] originalObservation = new int[li.size()];

		for (int m = 0; m < li.size(); m++) {
			originalObservation[m] = li.get(m);
		}

		int maxGeneCount = FindMaxArray.findMaxValueIntArray(originalObservation);
		int maxNodeSize = Math.max(100, maxGeneCount);

		if (maxNodeSize > 100) {
			for (Node n : queue) {
				n.setmaxNodeSize(maxNodeSize);
			}
		}
		ArrayList<Node> arln = new ArrayList<Node>(queue);
		return arln;
	}

	public static ArrayList<Node> setMaxNodeSize(Node root, int maxNodeSize) {

		Queue<Node> queue = root.postOrder();

		for (Node n : queue) {n.setmaxNodeSize(maxNodeSize);}
		ArrayList<Node> arln = new ArrayList<Node>(queue);
		return arln;
	}
}
