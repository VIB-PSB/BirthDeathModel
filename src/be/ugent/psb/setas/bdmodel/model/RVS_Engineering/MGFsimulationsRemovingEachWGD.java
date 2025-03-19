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
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;

public class MGFsimulationsRemovingEachWGD {

	private ArrayList<String> gfIDs;
	ProbCalculator probCalc;

	public ArrayList<String> getGfIDs() {return gfIDs;}

	public MGFsimulationsRemovingEachWGD(ArrayList<String> gfIDs,
			ProbCalculator probCache) {
		this.gfIDs = gfIDs;
		this.probCalc = probCache;
	}

	public MGFsimulationsRemovingEachWGD(ProbCalculator probCache) {
		this.gfIDs = new ArrayList<String>();
		this.probCalc = probCache;
	}

	public List<List<Double>> readMGFfile(String mgfFile) {

		FileReader fin = null;
		try {
			fin = new FileReader(mgfFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);
		// sc.nextLine(); /* The first line is a header: */

		List<List<Double>> MGFtable = new ArrayList<List<Double>>();
		gfIDs = new ArrayList<String>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();
			String[] chunks = line.split("\t");
			List<Double> nums = new ArrayList<Double>();
			gfIDs.add(chunks[0]);

			// i=0 GF-IDS ,i=1 rootSize, i=2 lambda, i=3 loglk, i=4 pv, the rest
			// gene counts
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

		// int[] mgfGeneCounts = new int[mgfRecord.size() - 4];
		int[] mgfGeneCounts = new int[mgfRecord.size() - 3];

		// for (int m = 4; m < mgfRecord.size(); m++) { /* rStar lamStar
		// loglkStar pvStar */
		for (int m = 3; m < mgfRecord.size(); m++) { /*
													 * rStar lamStar loglkStar
													 * pvStar
													 */
			mgfGeneCounts[m - 3] = mgfRecord.get(m).intValue();

		}

		return mgfGeneCounts;
	}

	public List<List<Integer>> getOnlyGFcounts(List<List<Double>> MGFtable) {

		List<List<Integer>> MGFintTable = new ArrayList<List<Integer>>();

		for (List<Double> mgfRecord : MGFtable) {

			List<Integer> ls = new ArrayList<Integer>();

			// for (int m = 4; m < mgfRecord.size(); m++) { // 0=r* 1=lam*
			// 2=loglk* 3=pv*
			for (int m = 3; m < mgfRecord.size(); m++) { // 0=r* 1=lam* 2=loglk*

				ls.add(mgfRecord.get(m).intValue());
			}
			MGFintTable.add(ls);
		}
		return MGFintTable;
	}

	public int getOnlyRootSize(List<Double> mgfRecord) {

		return mgfRecord.get(0).intValue();
	}

	public double calcLogLk_1MGF(Node root, List<Double> mgfRecord,
			ArrayList<Node> arln, ProbCalculator probCalc) {

		int rootStar = mgfRecord.get(0).intValue();
		double lambdaStar = mgfRecord.get(1);

		int[] mgfGeneCounts = new int[mgfRecord.size() - 3];

		for (int m = 3; m < mgfRecord.size(); m++) {
			mgfGeneCounts[m - 3] = mgfRecord.get(m).intValue();
		}
		root.setLeafValues(mgfGeneCounts);
		LikeLihood lkhd = new LikeLihood(lambdaStar, root.getmaxNodeSize() + 1,
				probCalc);
		double[] Lk_1MGF = lkhd.calcInternalLk(arln);
		double[] logLk_1MGF = MathOperations.giveLogArray(Lk_1MGF);

		return logLk_1MGF[rootStar];
	}

	public double calcCombinedLogLkAllMGF(Node root,
			List<List<Double>> MGFtable, ArrayList<Node> speciesTree) {

		double sumLoglk = 0;
		for (int i = 0; i < MGFtable.size(); i++) {
			List<Double> mgfRecord = MGFtable.get(i);
			double lk = calcLogLk_1MGF(root, mgfRecord, speciesTree,
					this.probCalc);
			sumLoglk += lk;
		}
		return sumLoglk;
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
//		double precisionLambda = 1e-5;
//		/** one digit more than the number of digits we require accurately */
//		double partitionSize = 0.1;
//		int defaultmaxNodeSize = 100;
//		// int numberOfObservations = 1000;
//
//		ProbCalculator probCalc = new ProbCalculator();
//		MGFsimulationsRemovingEachWGD mgfsim = new MGFsimulationsRemovingEachWGD(
//				probCalc);
//
//		String pathMGFfile = "/home/setas/Desktop/setas/Project1/RvsEng/Musa3TauMon2/gfAbsent/subSetMGFabsent";
////		List<List<Double>> MGFtable = mgfsim.readMGFfile(args[2]); // This table doesn't include the column with GF IDs: so later for reading GFcounts:only first 4 columns
//		
//		List<List<Double>> MGFtable = mgfsim.readMGFfile(pathMGFfile);
//		List<List<Integer>> gfCounts = mgfsim.getOnlyGFcounts(MGFtable);
//
////		int ignoreLineNumWGDFile = Integer.parseInt(args[3]);
//		int ignoreLineNumWGDFile =16;
//		/** fixed for each WGD in the .sh script */
////		int mgf_start = Integer.parseInt(args[4]);
//		int mgf_start =0;
//		/** PBS array ID we loop over */
////		 int mgf_end = Integer.parseInt(args[5]);
//		 
//		 int mgf_end =12;
//
////		int mgf_end = mgf_start + 800;
////
////		if (mgf_end > 9177) {
////			mgf_end = 9177;
////		}
////		int mgf =0;
//
//		
////		for(int ignoreLineNumWGDFile=0; ignoreLineNumWGDFile<=19; ignoreLineNumWGDFile++){
//		
//		Node root = SpeciesTreeParser.buildAndPartitionTree_ReverseEng(args[0],
//				args[1], partitionSize, defaultmaxNodeSize,
//				ignoreLineNumWGDFile);
//
//		 for (int mgf = mgf_start; mgf < mgf_end; mgf++) {
//		
//		 List<Double> mgfRecord = MGFtable.get(mgf);
//		 int rootSize = mgfRecord.get(0).intValue();
//		
//		 // System.out.print(mgfsim.getGfIDs().get(mgf) + "\t" + rootSize +
//		 // "\t" +mgfRecord.get(1)+"\t"+mgfRecord.get(2));
//		
//		 // for(int ignoreLineNumWGDFile=5; ignoreLineNumWGDFile<= 17;
//		 // ignoreLineNumWGDFile++){
//		
//		 ArrayList<Node> speciesTree = SpeciesTreeParser
//		 .setMaxNodeSizeAccToGF(root, gfCounts, mgf);
//		 SpeciesTreeParser.setLeavesValues(root, gfCounts, mgf);
//		
//		 // double lkAfterRemoval = mgfsim.calcLogLk_1MGF(root,
//		 // MGFtable.get(mgf), speciesTree);
//		
//		 // double lkAfterRemoval = mgfsim.calcCombinedLogLkAllMGF(root,
//		 // MGFtable,speciesTree);
//		 // System.out.print(ignoreLineNumWGDFile+ "\t"
//		 // + mgfRecord.get(0) + "\t" +
//		 // mgfRecord.get(1)+"\t"+lkAfterRemoval);
//		
//		 // System.out.println(lkAfterRemoval);
//		 // }
//		
//		 // System.out.println();
//		 // }
//		
//		 // for (int rootSize = 1; rootSize <= 20; rootSize++) {
//		
//		 ProbCalculator probCache = new ProbCalculator();
//		
//		 CuttingPlaneMethod cpm = new CuttingPlaneMethod(speciesTree,
//		 rootSize, stepSize, deltaLocalMoves, tolD, tolF,
//		 minInterval, maxInterval, minAllowed, maxAllowed,
//		 precisionLambda, probCache, root, gfCounts);
//		
//		 cpm.findOptimalLambda();
//		 
//		
//		 System.out.println(mgfsim.getGfIDs().get(mgf) + "\t" + rootSize
//		 + "\t" + mgfRecord.get(1) + "\t" + mgfRecord.get(2) + "\t"
//		 + cpm.getOptimalLambda() + "\t" + cpm.getFoptimalLambda());
//		
//		 // Pvalues pv = new Pvalues(root, speciesTree,
//		 // cpm.getOptimalLambda(), numberOfObservations, probCache);
//		 //
//		 // double pValue = pv.calculateConditionalPvalues(rootSize,
//		 // rootSize, cpm.getfOptimalLambda());
//		
//		  }
//
////		CommonFunctions cmmFunct = new CommonFunctions();
////		String idealProfilePath = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ComparisonsWithIdealprofiles/IdealProfiles_R1";
////		ArrayList<Integer> idealProfileR1 = cmmFunct.readColX_Int(idealProfilePath, 4); // Tau-Mon2_Musa3
////		
////		 ArrayList<Node> speciesTree = SpeciesTreeParser.setMaxNodeSize(root, 100);
////		 SpeciesTreeParser.setLeavesValues_one(root, idealProfileR1);
////		 
////	     double lkAfterRemoval = mgfsim.calcLogLk_1MGF(root, mgfRecord, speciesTree, probCalc);
////			 System.out.print(mgfsim.getGfIDs().get(mgf) + "\t" + mgfRecord.get(0) +"\t" +mgfRecord.get(1)+"\t"+mgfRecord.get(2)+"\t");
////	     
////	     System.out.println(lkAfterRemoval);
////		}
////	}
//
//}
}
