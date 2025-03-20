package utils.parsers;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import utils.bdmodel.Node;
import utils.bdmodel.PartitionBranches;

public class SpeciesTreeParser {
        
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
        
    public static Node buildInsertWGDsandPartitionTree_Negatives(String treeFile, String wgdFile, double partitionSize,
			int defaultmaxNodeSize) {

		NewickParser np = new NewickParser();
		Node root = np.buildTree(treeFile, defaultmaxNodeSize);
		root.getLeaves();

		WGMparser wgm = new WGMparser();
		List<List<String>> wgdList = wgm.readInputFile(wgdFile);
		root.insertWGMsToTheTree_Negatives(wgdList);

		if (partitionSize != 0) {
			PartitionBranches pb = new PartitionBranches(root, partitionSize);
			pb.insertAllVNSonAllBranches();
		}

		return root;
	}
	
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
        
    public static Node buildAndPartitionTree_ReverseEng_Negatives(String treeFile, String wgdFile, double partitionSize,
			int defaultmaxNodeSize, int includeLine) {

		NewickParser np = new NewickParser();
		Node root = np.buildTree(treeFile, defaultmaxNodeSize);
		root.getLeaves();

		WGMparser wgm = new WGMparser();
		List<List<String>> wgdList = wgm.readInputFile(wgdFile);
		root.insertWGMsToTheTree_ReverseEng_Negatives(wgdList, includeLine);

		if (partitionSize != 0) {
			PartitionBranches pb = new PartitionBranches(root, partitionSize);
			pb.insertAllVNSonAllBranches();
		}
		return root;
	}


	public static void setLeavesValues(Node root, List<List<Integer>> nonRepetetiveCounts, int gf) {

		List<Integer> li = nonRepetetiveCounts.get(gf);
		int[] originalObservation = new int[li.size()];

		for (int m = 0; m < li.size(); m++) {
			originalObservation[m] = li.get(m);
		}
		root.setLeafValues(originalObservation);
	}


	public static ArrayList<Node> setMaxNodeSize(Node root, int maxNodeSize) {

		Queue<Node> queue = root.postOrder();

		for (Node n : queue) {n.setmaxNodeSize(maxNodeSize);}
		ArrayList<Node> arln = new ArrayList<Node>(queue);
		return arln;
	}
}
