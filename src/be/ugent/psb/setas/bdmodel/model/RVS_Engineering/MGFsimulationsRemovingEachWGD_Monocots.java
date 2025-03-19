package be.ugent.psb.setas.bdmodel.model.RVS_Engineering;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ugent.psb.setas.bdmodel.model.CuttingPlaneMethod;
import be.ugent.psb.setas.bdmodel.model.LikeLihood;
import be.ugent.psb.setas.bdmodel.model.MathOperations;
import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.bdmodel.model.ProbCalculator;
import be.ugent.psb.setas.bdmodel.parsers.ReadGFcountsFile;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;

public class MGFsimulationsRemovingEachWGD_Monocots {

	private ArrayList<String> gfIDs;

	public ArrayList<String> getGfIDs() {
		return gfIDs;
	}

	public List<List<Double>> readMGFfile(String mgfFile) {

		FileReader fin = null;
		try {
			fin = new FileReader(mgfFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		/* The first line is a header: */
		Scanner sc = new Scanner(fin);
		sc.nextLine();

		List<List<Double>> MGFtable = new ArrayList<List<Double>>();
		gfIDs = new ArrayList<String>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();
			String[] chunks = line.split("\t");

			List<Double> nums = new ArrayList<Double>();
			gfIDs.add(chunks[0]);

			// i=0 GF-IDS ,i=1 rootSize, i=2 lambda, i=3 loglk, i=4 pv, the rest gene counts
			for (int i = 1; i < chunks.length; i++) { 	
				double parsed = Double.parseDouble(chunks[i]);
				nums.add(parsed);
			}

			MGFtable.add(nums);

		}
		sc.close();
		return MGFtable;

	}

	public int[] readMGFgeneCounts(List<Double> mgfRecord) {

		int[] mgfGeneCounts = new int[mgfRecord.size() - 4];

		for (int m = 4; m < mgfRecord.size(); m++) { // rStar lamStar loglkStar pvStar
			mgfGeneCounts[m - 4] = mgfRecord.get(m).intValue();

		}

		return mgfGeneCounts;
	}

	public List<List<Integer>> getOnlyGFcounts(List<List<Double>> MGFtable) {

		List<List<Integer>> MGFintTable = new ArrayList<List<Integer>>();

		for (List<Double> mgfRecord : MGFtable) {

			List<Integer> ls = new ArrayList<Integer>();

			for (int m = 4; m < mgfRecord.size(); m++) { //r* lam* loglk* pv* 

				ls.add(mgfRecord.get(m).intValue());
			}

			MGFintTable.add(ls);
		}

		return MGFintTable;
	}
	
	public double calcLogLk_1MGF(Node root,List<Double> mgfRecord, ArrayList<Node> arln) {

		int rootStar = mgfRecord.get(0).intValue();		
		double lambdaStar = mgfRecord.get(1);

		int[] mgfGeneCounts = new int[mgfRecord.size() - 4];
		
		for (int m = 4; m < mgfRecord.size(); m++) {
			mgfGeneCounts[m - 4] = mgfRecord.get(m).intValue();

		}

		root.setLeafValues(mgfGeneCounts);
		LikeLihood lkhd = new LikeLihood(lambdaStar,
				root.getmaxNodeSize() + 1);

		double[] Lk_1MGF = lkhd.calcInternalLk(arln);
		double [] logLk_1MGF = MathOperations.giveLogArray(Lk_1MGF);

		return logLk_1MGF[rootStar];
		
	}

//	public static void main(String[] args) {
//
//		double stepSize = 1e-4;// step calculating derivative
//		double deltaLocalMoves = 1e-1;
//		double tolD = 1e-3;
//		double tolF = 1e-4;
//		double minInterval = 1e-2;
//		double maxInterval = 10;
//		double minAllowed = 1e-2;
//		double maxAllowed = 10;
//		double precisionLambda = 1e-5; // one digit more than the number of digits we require accurately
//		double partitionSize = 0.1;
//		int defaultmaxNodeSize = 100;
////		int numberOfObservations = 1000;
//
//		MGFsimulationsRemovingEachWGD_Monocots mgfsim = new MGFsimulationsRemovingEachWGD_Monocots();
//		List<List<Double>> MGFtable = mgfsim.readMGFfile(args[2]); // Does not include the column with GF IDs: so later for reading GFcounts: remove only first 4 columns
//
//		List<List<Integer>> gfCounts = mgfsim.getOnlyGFcounts(MGFtable);
//
//		int ignoreLineNumWGDFile = Integer.parseInt(args[3]);
//		int mgf = Integer.parseInt(args[4]); //PBSarrayID or loop
//		
//		List<Double> mgfRecord = MGFtable.get(mgf);
//
//		  /** repeat the old data: */
//		  System.out.println(mgfsim.getGfIDs().get(mgf) + "\t"
//				+ mgfRecord.get(0) + "\t" + mgfRecord.get(1)+"\t"+mgfRecord.get(2)+"\t"+mgfRecord.get(3));
//			
//		   Node root = SpeciesTreeParser.buildInsertWGDsandPartitionTree(
//				args[0], args[1], partitionSize, defaultmaxNodeSize);
//		   ArrayList<Node> speciesTree = SpeciesTreeParser
//				.setMaxNodeSizeAccToGF(root, gfCounts, mgf);
//		   SpeciesTreeParser.setLeavesValues(root, gfCounts, mgf);
//		   
//		   double lkFixedRlamdaMonocots = mgfsim.calcLogLk_1MGF(root, MGFtable.get(mgf), speciesTree);
//			System.out.println("fixRLamStar\t"+ mgfRecord.get(0) + "\t" + mgfRecord.get(1)+"\t"+lkFixedRlamdaMonocots);
//		
//		    //build the tree again, add WGDs except one, and partition the branches (To refresh all previuos settings)
//			Node rootAfterRemoval1WGD = SpeciesTreeParser.buildAndPartitionTree_ReverseEng(
//					args[0], args[1], partitionSize, defaultmaxNodeSize,
//					ignoreLineNumWGDFile);
//			ArrayList<Node> speciesTreeAfterRemoval1WGD = SpeciesTreeParser
//					.setMaxNodeSizeAccToGF(rootAfterRemoval1WGD, gfCounts, mgf);
//			SpeciesTreeParser.setLeavesValues(rootAfterRemoval1WGD, gfCounts, mgf);
//			
//			
//			ProbCalculator probCache1 = new ProbCalculator();
//			int rStar = mgfRecord.get(0).intValue();
//			
//			CuttingPlaneMethod cpm1 = new CuttingPlaneMethod(speciesTree,
//					rStar, stepSize, deltaLocalMoves, tolD, tolF,
//					minInterval, maxInterval, minAllowed, maxAllowed,
//					precisionLambda, probCache1);
//			
//			cpm1.findOptimalLambda();
//			//fixing the rStar what is the lambdaStar in the tree of Monocots?
//			System.out.println("fixRCPMlam"
//					+ ""
//					+ ""+"\t"+rStar + "\t" + cpm1.getOptimalLambda()
//					+ "\t" + cpm1.getFoptimalLambda());
//			// we can here see if the lambdaStar(MGF) for Monocots is different than the one for all Angiosperms. (higher or lower?)
//			// we can then compare this to the optimal lambda of the tree without 1WGD (which might be obtained for a different rootSize)
//			
//			
//			/** after removal of 1 WGD */
//			
//			double lkAfterRemoval = mgfsim.calcLogLk_1MGF(rootAfterRemoval1WGD, MGFtable.get(mgf), speciesTreeAfterRemoval1WGD);
//			System.out.print(ignoreLineNumWGDFile+ "\t"
//					+ mgfRecord.get(0) + "\t" + mgfRecord.get(1)+"\t"+lkAfterRemoval);
//			
//			if(lkAfterRemoval < lkFixedRlamdaMonocots){System.out.print("\t"+ "lower");}
//			
//			else{System.out.print("\t"+ "higher");}
//			
//			System.out.println();
//			
//			for (int rootSize = 1; rootSize <= 20; rootSize++) {
//				
//				// To refresh all the settings for the new calculations for the new root size
//				speciesTreeAfterRemoval1WGD = SpeciesTreeParser
//						.setMaxNodeSizeAccToGF(rootAfterRemoval1WGD, gfCounts, mgf);	
//				SpeciesTreeParser.setLeavesValues(rootAfterRemoval1WGD, gfCounts, mgf);
//
//				ProbCalculator probCache = new ProbCalculator();
//
//				CuttingPlaneMethod cpm = new CuttingPlaneMethod(speciesTreeAfterRemoval1WGD,
//						rootSize, stepSize, deltaLocalMoves, tolD, tolF,
//						minInterval, maxInterval, minAllowed, maxAllowed,
//						precisionLambda, probCache);
//
//				cpm.findOptimalLambda();
//				System.out.print(rootSize + "\t" + cpm.getOptimalLambda()
//						+ "\t" + cpm.getFoptimalLambda());
////				Pvalues pv = new Pvalues(root, speciesTree,
////						cpm.getOptimalLambda(), numberOfObservations, probCache);
////
////				double pValue = pv.calculateConditionalPvalues(rootSize,
////						rootSize, cpm.getfOptimalLambda());
////				System.out.print("\t" + pValue);
//				System.out.println();
//		}
//	}

}

