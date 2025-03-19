package be.ugent.psb.setas.bdmodel.model.RVS_Engineering;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.bdmodel.model.LikeLihood;
import be.ugent.psb.setas.bdmodel.model.MathOperations;
import be.ugent.psb.setas.bdmodel.model.Node;

public class inferOneLambdaAllMGFsForTree {
	
	public double[] calcMatrix(ArrayList<Node> arln, int columns,
			double stepSize) {
        double [] lkRootSize1 = new double[1000];
		for (int l = 0; l < 1000; l++) {
			double lambda = l * stepSize;
			LikeLihood lk = new LikeLihood(lambda, 100);
			lkRootSize1[l] = lk.calcInternalLk(arln)[1];
			lkRootSize1[l] = lk.calcInternalLk(arln)[2];
			lkRootSize1[l] = lk.calcInternalLk(arln)[3];	
		}
		return lkRootSize1;
	}
	
	public double [] calcLogLk_1MGF(Node root,List<Double> mgfRecord, ArrayList<Node> arln, double lambda, int rootSize) {

		int[] mgfGeneCounts = new int[mgfRecord.size() - 4];
		
		for (int m = 4; m < mgfRecord.size(); m++) {
			mgfGeneCounts[m - 4] = mgfRecord.get(m).intValue();

		}

		root.setLeafValues(mgfGeneCounts);
		LikeLihood lkhd = new LikeLihood(lambda,root.getmaxNodeSize() + 1);

		double[] Lk_1MGF = lkhd.calcInternalLk(arln);
		double [] logLk_1MGF = MathOperations.giveLogArray(Lk_1MGF);

		return logLk_1MGF;		
	}
	
//	public static void main(String [] args){
//		
//		inferOneLambdaAllMGFsForTree inferlam = new inferOneLambdaAllMGFsForTree();
//		double partitionSize =0.1;
//		int  defaultmaxNodeSize=100;
//		
//		MGFsimulationsRemovingEachWGD mgfsim = new MGFsimulationsRemovingEachWGD();
//		List<List<Double>> MGFtable = mgfsim.readMGFfile(args[2]); // Does not include the column with GF IDs: so later for reading GFcounts: remove only first 4 columns
//
//		List<List<Integer>> gfCounts = mgfsim.getOnlyGFcounts(MGFtable);
//
//		int ignoreLineNumWGDFile = Integer.parseInt(args[3]);
//		int mgf = Integer.parseInt(args[4]); //PBSarrayID or loop
//		
//		List<Double> mgfRecord = MGFtable.get(mgf);
//
//		System.out.println(mgfsim.getGfIDs().get(mgf) + "\t"
//				+ mgfRecord.get(0) + "\t" + mgfRecord.get(1)+"\t"+mgfRecord.get(2)+"\t"+mgfRecord.get(3)+"\t");
//		
//		
//		Node root = SpeciesTreeParser.buildAndPartitionTree_ReverseEng(
//				args[0], args[1], partitionSize, defaultmaxNodeSize,
//				ignoreLineNumWGDFile);
//
//		ArrayList<Node> speciesTree = SpeciesTreeParser
//				.setMaxNodeSizeAccToGF(root, gfCounts, mgf);
//
//		SpeciesTreeParser.setLeavsValues(root, gfCounts, mgf);
//		
//		int rootSize =1;
//		double lambda = 0.1;
//		double [] lkAfterRemoval = inferlam.calcLogLk_1MGF(root, mgfRecord, speciesTree, lambda, rootSize);
//		
//		
//		System.out.print(ignoreLineNumWGDFile+ "\t"
//				+ mgfRecord.get(0) + "\t" + mgfRecord.get(1)+"\t"+lkAfterRemoval);
//	}
	
	
	
	
	
	
}
