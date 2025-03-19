package be.ugent.psb.setas.bdmodel.model.RVS_Engineering;

//import be.ugent.psb.setas.bdmodel.model.GenerateObservations;
//import be.ugent.psb.setas.bdmodel.model.Node;
//import be.ugent.psb.setas.bdmodel.model.ProbCalculator;
//import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;
//import be.ugent.psb.setas.independent_parsers.CommonFunctions;
//
public class BruteForceRecheckLRTs {
//
//	public int[] findObservationBasedOnGFid(String gfcountfile, String gfID, int numberOfLeaves) {
//
//		int[] observation = new int[numberOfLeaves];
//		CommonFunctions comFunc = new CommonFunctions();
//		List<List<String>> map = comFunc.readMapFile(gfcountfile);
//
//		for (int i = 0; i < map.size(); i++) {
//			if (map.get(i).get(0).equals(gfID)) {
//				for (int j = 0; j < observation.length; j++) {
//					observation[j] = Integer.parseInt(map.get(i).get(j + 3));
//				}
//				return observation;
//			}
//		}
//
//		return observation;
//	}
//
//	public int calculateSumStatisticsOfBruteForceResults(String resultsFile, int gfNumberInList) {
//
//		CommonFunctions cmf = new CommonFunctions();
//		List<List<String>> allResults = cmf.readMapFile(resultsFile);
//		List<String> gfRecord = allResults.get(gfNumberInList);
//
//		double[] lrts = new double[1000];
//		double lrtOriginal = Double.parseDouble(gfRecord.get(3));
//
//		int sum = 0;
//
//		for (int i = 4; i < gfRecord.size(); i++) {
//
//			lrts[i] = Double.parseDouble(gfRecord.get(i));
//
//			if (lrts[i] > lrtOriginal) {
//
//				sum += 1;
//			}
//
//		}
//
//		return sum;
//
//	}
//
//	public static void main(String[] args) {
//		/**
//		 * if H0= full-tree and H1= Rm-WGDi, all r* and lambda* in all loglk //
//		 * values must come from H0, so: // read the GFid r* lambda* from file
//		 * sortedDefault2_rmWGD0, fix r // and lambda to r* and lambda*, remove
//		 * WGD0 and recalculate // loglk(rmWGD0|r* and lam*) //* print GF-id r*
//		 * lam* loglk* loglk*(tree-WGD0| r* and lam*) // // we will have 20
//		 * lines like this for each WGD, so we better print // // GF-id r* lam*
//		 * 2(loglk*(tree-WGD0| r* and lam*) - loglk* ) // =LRT-original*
//		 **/
//
//		int numberWGDs = 20;
//		double partitionSizeBranch = 0.1;
//		int defaultmaxNodeSize = 100;
//
//		int lengthMCMC = 1000;
//		int numOfObservations = 1000;
//		//
//		CalculateCorrectLRToriginalFiles calcLrt = new CalculateCorrectLRToriginalFiles();
//		CommonFunctions cmf = new CommonFunctions();
//		ProbCalculator probCalc = new ProbCalculator();
//		//
//		// String combinedOutput = args[6];
//		// String gfCountsFile = args[7];
//		String combinedOutput = "/home/setas/Desktop/setas/Project1/Results/CompareRankings/newCombinedOutput/combinedOutputFiles_rawSorted/SortedLambda/37spe/comOut_37spe_Tau_Mon2_Musa3_SortedLambda";
//		String gfCountsFile = "/home/setas/Desktop/setas/Project1/Simulations/37speMGCF5-CPMpval-Tier1/37spe-MGCF5-9178core-OrderNewickTree-noHeader";
//
//		Node rootOriginal = SpeciesTreeParser.buildInsertWGDsandPartitionTree(args[0], args[1], partitionSizeBranch,
//				defaultmaxNodeSize);
//		GenerateObservations go = new GenerateObservations(0, 100, false, probCalc, lengthMCMC);
//
//		List<List<String>> map = cmf.readMapFile(combinedOutput);
//		for (int i = 0; i < map.size(); i++) {
//
//			List<String> gfrecord = map.get(i);
//			String gfID = gfrecord.get(0);
//
//			int rootSizeStar = Integer.parseInt(gfrecord.get(1));
//			double lambdaStar = Double.parseDouble(gfrecord.get(2));
//			double loglkStar = Double.parseDouble(gfrecord.get(3));
//
//			int sum = 0;
//
//			System.out.print(gfID + "\t" + rootSizeStar + "\t" + lambdaStar);
//
//			// lambdaPartitionSize =0.1, not needed here
//			GenerateBackgroundDistLRT genBg = new GenerateBackgroundDistLRT(0.1, rootSizeStar, probCalc);
//			ArrayList<Node> TreesWithoutWGDs = genBg.rootRmWGDs(args[0], args[1], args[2], args[3], args[4], args[5],
//					numberWGDs, partitionSizeBranch, defaultmaxNodeSize);
//
//			int[] observation_original = calcLrt.findObservationBasedOnGFid(gfCountsFile, gfID, 37);
//
//			// for (int rmWGD = 0; rmWGD < 20; rmWGD++) {
//			int rmWGD = 1;
//
//			Node rootRmWGDi = TreesWithoutWGDs.get(rmWGD);
//			double loglkRmWGDi = genBg.calcLoglk(rootRmWGDi, lambdaStar, observation_original);
//
//			double lrt_original = 2 * (loglkRmWGDi - loglkStar);
//			System.out.print("\t" + lrt_original);
//			// }
//
//			for (int obs = 0; obs < numOfObservations; obs++) {
//
//				int[] observation_H0 = go.generateObservation(rootOriginal, rootSizeStar, lambdaStar);
//				double loglkRmWGDi_H0 = genBg.calcLoglk(rootRmWGDi, lambdaStar, observation_H0);
//				double loglkOriginal_H0 = genBg.calcLoglk(rootOriginal, lambdaStar, observation_H0);
//
//				double lrt_simulated = 2 * (loglkRmWGDi_H0 - loglkOriginal_H0);
//
//				// System.out.print("\t"+lrt_simulated);
//
//				if (lrt_simulated > lrt_original) {
//					sum += 1;
//				}
//			}
//
//			System.out.print("\t" + sum);
//			System.out.print("\n");
//		}
//
//		/** if H0 = rmWGDi , H1=Full-tree **/
//		// String RVS_rm20WGDs =
//		// "/home/setas/Desktop/setas/Project1/Simulations/LRTbackground/sorted_TauMon2Musa3_RmWGDs.txt";
//		// String RVS_rm20WGDs = args[8];
//		// List<List<String>> map = cmf.readMapFile(RVS_rm20WGDs);
//		//
//		// int start=0;
//		// int end = 4500;
//		/**
//		 * because of time limits in tier2. we have to repeat for up to end=
//		 * 9178
//		 **/
//
//		/** For incompelte tier2 output : new jar file! **/
//		// String output_tier2_incomplete = args[9];
//		// String output_tier2_incomplete
//		// ="/home/setas/Desktop/setas/Project1/RvsEng/Musa3TauMon2/H0isRmWGD/output_rmWGDi_original_merged";
//		// ArrayList<String> gfIDsIntermediateTier2Output =
//		// cmf.readColX_String(output_tier2_incomplete, 0);
//
//		// for (int i = start; i <= end; i++) {
//		//
//		// List<String> gfrecord = map.get(i);
//		// String gfID = gfrecord.get(0);
//		//
//		// if(!gfIDsIntermediateTier2Output.contains(gfID)){
//		// int rootSizeStar = Integer.parseInt(gfrecord.get(1));
//		//
//		// GenerateBackgroundDistLRT genBg = new
//		// GenerateBackgroundDistLRT(0.1,rootSizeStar,probCalc);
//		//
//		// int[] observation = calcLrt.findObservationBasedOnGFid(gfCountsFile,
//		// gfID, 37);
//		//
//		// double [] lambdas_rmWGDs = new double[20]; //this is optimzied after
//		// removal of the wgd
//		// double [] loglks_rmWGDs = new double[20];
//		//
//		// for(int m=4; m<44 ;m+=2){
//		// lambdas_rmWGDs[(m-4)/2]= Double.parseDouble(gfrecord.get(m));
//		// loglks_rmWGDs[(m-4)/2]=Double.parseDouble(gfrecord.get(m+1));
//		// }
//		//
//		// Node rootFullTree = SpeciesTreeParser.buildAndPartitionTree(args[0],
//		// args[1], partitionSizeBranch, defaultmaxNodeSize);
//		// rootFullTree.setLeafValues(observation);
//		//
//		// System.out.print(gfID+"\t"+rootSizeStar); // we skip printing lambda
//		// this time, because it is 20 lam per family, so we have to adjust the
//		// python script to skip one column
//		//
//		// for(int k=0; k<numberWGDs; k++){
//		//
//		// double loglkFullTree = genBg.calcLoglk(rootFullTree,
//		// lambdas_rmWGDs[k], observation);
//		//
//		// double lrt = 2*(loglkFullTree-loglks_rmWGDs[k]);
//		// System.out.print("\t"+lrt);
//		// }
//		//
//		// System.out.print("\n");
//		//
//		// }
//		// }
//
//		// BruteForce teste results: to see if it is different than tier2
//		// results for categorized lambda
//
//		String resultFile = "/home/setas/Desktop/lrts_BruteForce";
//		BruteForceRecheckLRTs brtforce = new BruteForceRecheckLRTs();
//
//		CommonFunctions cmf2 = new CommonFunctions();
//		List<List<String>> allResults = cmf2.readMapFile(resultFile);
//
//		for (int i = 0; i < 475; i++) {
//
//			List<String> gfRecord = allResults.get(i);
//
//			System.out.print(gfRecord.get(0) + "\t");
//
//			int sum = brtforce.calculateSumStatisticsOfBruteForceResults(resultFile, i);
//
//			System.out.print("\t" + sum);
//		}
//
//		System.out.print("\n");
//
//	}
}
