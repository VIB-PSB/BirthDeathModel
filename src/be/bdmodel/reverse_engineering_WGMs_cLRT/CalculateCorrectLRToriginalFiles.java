package be.ugent.psb.setas.bdmodel.model.RVS_Engineering;

import java.util.List;
import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class CalculateCorrectLRToriginalFiles {

	public int[] findObservationBasedOnGFid(String gfcountfile, String gfID, int numberOfLeaves) {

		int[] observation = new int[numberOfLeaves];
		CommonFunctions comFunc = new CommonFunctions();
		List<List<String>> map = comFunc.readMapFile(gfcountfile);

		for (int i = 0; i < map.size(); i++) {
			if (map.get(i).get(0).equals(gfID)) {
				for (int j = 0; j < observation.length; j++) {
					observation[j] = Integer.parseInt(map.get(i).get(j + 3));
				}
				return observation;
			}
		}

		return observation;
	}
	

//	public static void main(String[] args) {
////		/** if H0= full-tree and H1= Rm-WGDi, all r* and lambda* in all loglk values must come from H0, so:
////		// read the GFid r* lambda* from file sortedDefault2_rmWGD0, fix r and
////		// lambda to r* and lambda*, remove WGD0 and recalculate loglk(rmWGD0|r* and lam*)
////		// print GF-id r* lam* loglk* loglk*(tree-WGD0| r* and lam*)
////		// we will have 20 lines like this for each WGD, so we better print
////		// GF-id r* lam* 2(loglk*(tree-WGD0| r* and lam*) - loglk* ) =LRT-original*
////		**/
////		
//		int numberWGDs = 20;
//		double partitionSizeBranch = 0.1;
//		int defaultmaxNodeSize = 100;
////
//		CalculateCorrectLRToriginalFiles calcLrt = new CalculateCorrectLRToriginalFiles();
//		CommonFunctions cmf = new CommonFunctions();
//		ProbCalculator probCalc = new ProbCalculator();
////
////		String combinedOutput = args[6];
////		String gfCountsFile = args[7];
//		String combinedOutput = "/home/setas/Desktop/setas/Project1/Results/CompareRankings/newCombinedOutput/combinedOutputFiles_rawSorted/SortedLambda/37spe/comOut_37spe_Tau_Mon2_Musa3_SortedLambda";
//		String gfCountsFile = "/home/setas/Desktop/setas/Project1/Simulations/37speMGCF5-CPMpval-Tier1/37spe-MGCF5-9178core-OrderNewickTree-noHeader";
//
////		List<List<String>> map = cmf.readMapFile(combinedOutput);	
////		for (int i = 0; i < map.size(); i++) {
////
////			List<String> gfrecord = map.get(i);
////			String gfID = gfrecord.get(0);
////
////			int rootSizeStar = Integer.parseInt(gfrecord.get(1));
////			double lambdaStar = Double.parseDouble(gfrecord.get(2));
////			double loglkStar = Double.parseDouble(gfrecord.get(3));
////			
////			//lambdaPartitionSize =0.1, not needed here
////			GenerateBackgroundDistLRT genBg = new GenerateBackgroundDistLRT(0.1,rootSizeStar,probCalc);
////			ArrayList<Node> TreesWithoutWGDs = genBg.rootRmWGDs(args[0], args[1], args[2], args[3], args[4], args[5],
////					numberWGDs, partitionSizeBranch, defaultmaxNodeSize);
////
////			int[] observation = calcLrt.findObservationBasedOnGFid(gfCountsFile, gfID, 37);			
////			System.out.print(gfID+"\t"+rootSizeStar+"\t"+lambdaStar);
////			
////			for(int rmWGD=0; rmWGD<20; rmWGD++){
////				
////			Node rootRmWGDi = TreesWithoutWGDs.get(rmWGD);
////			double loglkRmWGD0 = genBg.calcLoglk(rootRmWGDi, lambdaStar, observation);
////			
////			double lrt = 2*(loglkRmWGD0-loglkStar);
////			System.out.print("\t"+lrt);}			
////			
////			System.out.print("\n");
////		}
//				
//		/** if H0 = rmWGDi , H1=Full-tree**/
//		String RVS_rm20WGDs = "/home/setas/Desktop/setas/Project1/Simulations/LRTbackground/sorted_TauMon2Musa3_RmWGDs.txt";
////		String RVS_rm20WGDs = args[8];	
//		
//		List<List<String>> map = cmf.readMapFile(RVS_rm20WGDs);
//		
//		int start=5206; //5206
//		int end = 9178;	
//		/**because of time limits in tier2. we have to repeat for up to end= 9178**/
//		
//		/** For incompelte tier2 output : new jar file!**/
////		String output_tier2_incomplete = args[9];
////		String output_tier2_incomplete ="/home/setas/Desktop/setas/Project1/RvsEng/Musa3TauMon2/H0isRmWGD/output_rmWGDi_original_merged";
////		ArrayList<String> gfIDsIntermediateTier2Output = cmf.readColX_String(output_tier2_incomplete, 0); 
//
////		for (int i = start; i <= end; i++) {
////			
////		List<String> gfrecord = map.get(i);
////		String gfID = gfrecord.get(0);
////		
//////		if(!gfIDsIntermediateTier2Output.contains(gfID)){
////		int rootSizeStar = Integer.parseInt(gfrecord.get(1));
////		
////		GenerateBackgroundDistLRT genBg = new GenerateBackgroundDistLRT(0.1,rootSizeStar,probCalc);
////		
////		int[] observation = calcLrt.findObservationBasedOnGFid(gfCountsFile, gfID, 37);
////		
////		double [] lambdas_rmWGDs = new double[20]; //this is optimzied after removal of the wgd
////		double [] loglks_rmWGDs = new double[20];
////		
////		for(int m=4; m<44 ;m+=2){
////		lambdas_rmWGDs[(m-4)/2]= Double.parseDouble(gfrecord.get(m));
////		loglks_rmWGDs[(m-4)/2]=Double.parseDouble(gfrecord.get(m+1));	
////		}
////		
////		Node rootFullTree = SpeciesTreeParser.buildInsertWGDsandPartitionTree(args[0], args[1], partitionSizeBranch, defaultmaxNodeSize);
////		rootFullTree.setLeafValues(observation);
////		
////		System.out.print(gfID+"\t"+rootSizeStar); // we skip printing lambda this time, because it is 20 lam per family, so we have to adjust the python script to skip one column and read 20 diff lambdas
////				
////		for(int k=0; k<numberWGDs; k++){			
////		
////		double loglkFullTree = genBg.calcLoglk(rootFullTree, lambdas_rmWGDs[k], observation);
////		
////		double lrt = 2*(loglkFullTree-loglks_rmWGDs[k]);
////		System.out.print("\t"+lrt);
////		}
////		
////		System.out.print("\n");
////		
//////	}
////		}
//		
//		
//		
//	}
}
