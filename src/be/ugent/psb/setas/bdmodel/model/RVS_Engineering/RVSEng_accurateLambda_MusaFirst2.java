package be.ugent.psb.setas.bdmodel.model.RVS_Engineering;

import java.util.ArrayList;

import be.ugent.psb.setas.bdmodel.model.CuttingPlaneMethod;
import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.bdmodel.model.ProbCalculator;

public class RVSEng_accurateLambda_MusaFirst2 {
	
	final double stepSize = 1e-4;// step calculating derivative
	final double deltaLocalMoves = 1e-1;
	final double tolD = 1e-3;
	final double tolF = 1e-4;
	final double minInterval = 1e-2;
	final double maxInterval = 9.99;
	final double minAllowed = 1e-2;
	final double maxAllowed = 9.99;
	final double precisionLambda = 1e-5;
	
	private int maxRootNodeSize;
	
	public RVSEng_accurateLambda_MusaFirst2(int rootSize) {
		this.maxRootNodeSize = rootSize;
	}
		
//	public double[] optimize_r_lambda(ArrayList<Node> speTree, Node root, ProbCalculator probCalc) {
//
//		double[] rStar_lambdaStar_LoglkStar = new double[3];
//		
//		double[][] r_lambda_lkMatrix = new double[10][99999];
//			
//			int rowNumber =0;
//			for (int lam_col = 1; lam_col < 100000; lam_col++) {
//				double testLambda = lam_col * (0.0001);
//				LikeLihood lk = new LikeLihood(testLambda, speTree.get(0).getmaxNodeSize() + 1, probCalc);
//				r_lambda_lkMatrix[rowNumber]  = lk.calcInternalLk(speTree);
//			}
//		return rStar_lambdaStar_LoglkStar;
//	}
	
	public double[] cpm_optimize_r_lambda(ArrayList<Node> speTree, Node root, ProbCalculator probCalc) {
		
		double[] rStar_lambdaStar_LoglkStar = new double [3];
		
		for(int testRootSizes = 1; testRootSizes <= 10; testRootSizes++) {
			
			root.setValue(testRootSizes);
			
			double[] tmp= cpm_fixedRootSize (speTree, testRootSizes, root, probCalc);
			
			if(testRootSizes==1) {
				rStar_lambdaStar_LoglkStar = tmp;}
			
			else { //find the rootSize with maximum loglk
				if(tmp[2] > rStar_lambdaStar_LoglkStar[2]) {
					rStar_lambdaStar_LoglkStar = tmp;			
				}
			}
		}
		
		return rStar_lambdaStar_LoglkStar;
	}
	
	public double[] cpm_fixedRootSize (ArrayList<Node> speTree, int rootSize, Node root, ProbCalculator probCalc){
		
		double[] rStar_lamStar_LoglkStar = new double [3];
		rStar_lamStar_LoglkStar[0] = rootSize; 
		
		 CuttingPlaneMethod cpm = new CuttingPlaneMethod(speTree,
		 rootSize, stepSize, deltaLocalMoves, tolD, tolF,
		 minInterval, maxInterval, minAllowed, maxAllowed,
		 precisionLambda, probCalc, root);	 
		 cpm.findOptimalLambda();	
		 rStar_lamStar_LoglkStar[1] = cpm.getOptimalLambda();
		 rStar_lamStar_LoglkStar[2] = cpm.getFoptimalLambda();
		 return rStar_lamStar_LoglkStar;
	}
	

//	public static void main(String[] args) {
//
//		int defaultmaxNodeSize = 100;
//		double partitionSize = 0.1;
//		int lengthMCMC = 1000;
//		int numOfObservations = 1000;
//		int numberWGDs = 20;
//		
//		CommonFunctions cmf = new CommonFunctions();
//		ProbCalculator probCalc = new ProbCalculator();
//		GenerateBackgroundDistLRT genBgLRT = new GenerateBackgroundDistLRT(probCalc);
//		ReadGFcountsFile rgf = new ReadGFcountsFile();
//		CalculateCorrectLRToriginalFiles calcLrt = new CalculateCorrectLRToriginalFiles();
//		RVSEng_accurateLambda_MusaFirst2 rvsAcc = new RVSEng_accurateLambda_MusaFirst2();
//		
////		String treeFile = args[0];
////		String wgdFile = args[1];
////		String wgdFile_rmWGD8 = args[2];
////		String wgdFile_rmWGD9=args[3];
////		String wgdFile_rmWGD12=args[4];
////		String wgdFile_rmWGD16=args[5];
////		String combinedOutputFile = args[6];
////		String gfCountsFile = args[7];
//
//
//		String treeFile = "/BirthDeathModel/data/bdfiles/rvs_new/37spe-MGCF5.txt";
//		String wgdFile = "/BirthDeathModel/data/bdfiles/rvs_new/WGDs-MusaFirst2Successive";
//		String wgdFile_rmWGD8 = "/BirthDeathModel/data/bdfiles/rvs_new/WGDs-MusaFirst2Successive_rm8";
//		String wgdFile_rmWGD9 = "/BirthDeathModel/data/bdfiles/rvs_new/WGDs-MusaFirst2Successive_rm9";
//		String wgdFile_rmWGD12 = "/BirthDeathModel/data/bdfiles/rvs_new/WGDs-MusaFirst2Successive_rm12";
//		String wgdFile_rmWGD16 = "/BirthDeathModel/data/bdfiles/rvs_new/WGDs-MusaFirst2Successive_rm16";
//		String combinedOutputFile = "/BirthDeathModel/data/bdfiles/rvs_new/combinedOutput_MusaFirst2close";
//		String gfCountsFile = "/BirthDeathModel/data/bdfiles/rvs_new/37spe-MGCF5-9178core-OrderNewickTree-noHeader";
//		
//		List<List<Integer>> gfCounts = rgf.read_all(gfCountsFile);
//		
//		Node rootOriginal = SpeciesTreeParser.buildInsertWGDsandPartitionTree(treeFile, wgdFile, partitionSize,
//				defaultmaxNodeSize);
//
//
//		ArrayList<Node> rootRmWGDs = genBgLRT.rootRmWGDs(treeFile, wgdFile, wgdFile_rmWGD8, wgdFile_rmWGD9, wgdFile_rmWGD12, wgdFile_rmWGD16,
//				numberWGDs, partitionSize, defaultmaxNodeSize);
//		
//		Queue<Node> nodes = rootRmWGDs.get(8).postOrder();
//		
//		/** read the gene family ID, rootSizeStar, lambdaStar, loglkStar under the null model*/
//	
//		List<List<String>> combinedOutput = cmf.readMapFile(combinedOutputFile);
//		
//		int gfNumber = 0;
////		int gfNumber =Integer.parseInt(args[8]);
//
//		List<String> myGFrecord = combinedOutput.get(gfNumber);
//		String gfID = myGFrecord.get(0);
//		int rootSizeStar = Integer.parseInt(myGFrecord.get(1));
//		double lambdaStar_fullTree = Double.parseDouble(myGFrecord.get(2));
//		double loglkStar_fullTree = Double.parseDouble(myGFrecord.get(3));
//
//		System.out.print(gfID + "\t" + rootSizeStar + "\t" + lambdaStar_fullTree + "\t" + loglkStar_fullTree+"\t");
//
//		int[] originalGFCounts = calcLrt.findObservationBasedOnGFid(gfCountsFile, gfID, 37);
//
//		int numberOfWGD = 0;
////		for(int numberOfWGD = 0; numberOfWGD < 20; numberOfWGD++) {
//
//			Node rootRmWGD = rootRmWGDs.get(numberOfWGD);
//
//			rootRmWGD.setLeafValues(originalGFCounts);
//
//			ArrayList<Node> speTree_rmWGDi = SpeciesTreeParser.setMaxNodeSize(rootRmWGD, defaultmaxNodeSize);
//
//			double[] lambdaStar_LoglkStar_rm = rvsAcc.CPM_optimal_R_lambda(speTree_rmWGDi, rootRmWGD, probCalc);
//
//			for (double d : lambdaStar_LoglkStar_rm) {
//				System.out.print("\t" + d);
//			}
//
//			System.out.print("\n");
////		}
//		
//		
//
//// For this part, we need to get the new optimal lambdas for rm_WGDs first.
////		double[] loglkRmWGDs = new double[numberWGDs];
////		double[] lrtRmWGDs = new double[numberWGDs];
////		GenerateObservations go = new GenerateObservations(0, 100, false, probCalc, lengthMCMC);
////
////		for (int i = 0; i < numOfObservations; i++) {
////			// for(int i=0; i<numberWGDs;i++){
////
////			int[] observation_H0 = go.generateObservation(rootRmWGDs.get(i), rootSizeStar, lambdaStar);
////
////			loglkRmWGDs[i] = genBgLRT.calcLoglk(rootRmWGDs.get(i), lambdaStar, observation_H0);
////
////			double loglkFullTree = genBgLRT.calcLoglk(rootOriginal, lambdaStar, observation_H0);
////
////			lrtRmWGDs[i] = 2 * (loglkFullTree - loglkRmWGDs[i]);
////
////			System.out.print(lrtRmWGDs[i] + "\t");
////
////		}
//
//	}
}
