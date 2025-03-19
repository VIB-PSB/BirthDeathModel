package be.ugent.psb.setas.independent_parsers.Sort;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class SortAllGFAllGenesInOrderLam {

	public static void main(String[] args) {
		CommonFunctions cmmFunct = new CommonFunctions();

		String pathGFIDsInOrder = "/home/setas/Desktop/group_esb/setas_stmae_locar/Reviews/AllGFs_InOrderLam_MusaFirst2Close_PartitionSC";
		ArrayList<String> gfIDs_SC = cmmFunct.readColX_String(pathGFIDsInOrder, 0);
		ArrayList<Integer> groupID = cmmFunct.readColX_Int(pathGFIDsInOrder, 1);
		
//		String topGFidsFile= "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/MusaFirst2WGDs/basedOnCombinedRanking/TopBottom/Top1000GFids_combRank";
//		String bottomGFidsFile="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/MusaFirst2WGDs/basedOnCombinedRanking/TopBottom/Bottom1000GFids_combRank";	
//		List<String> topGFs= cmmFunct.read1ColFile_String(topGFidsFile);
//		List<String> bottomGFs=cmmFunct.read1ColFile_String(bottomGFidsFile);

//		String pathGFIDsGeneIDs = "/home/setas/Desktop/KnKsT_AllGenePairs_Allspe_MusaFisrt2Close-CombRank";
//		String pathGFIDsGeneIDs ="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/Exp_Func_newTopBottom1000/genefamily.angiosperm.core.Atha.codeml.Coex.TopBottomMiddle";
//		String pathGFIDsGeneIDs="/home/setas/Desktop/TopBottomMiddle_CombRank_MusaFirst2Close_ExpFunctPlots";
//		List<List<String>> GFGeneMap = cmmFunct.readMapFile(pathGFIDsGeneIDs);

		boolean found = false;

		// for (int i = 0; i < gfIDs.size(); i++) {
		//
		// String gfID = gfIDs.get(i);
		//
		// for (List<String> ls : GFGeneMap) {
		//
		// if (gfID.equals(ls.get(0))) {
		//
		// for (String s : ls) {
		// System.out.print(s + "\t");
		// }
		//
		// System.out.print(groupID.get(i) + "\n");
		// }
		// }
		// }

//		for (int i = 1; i < GFGeneMap.size(); i++) { //expFunc source file has a header
//
//			List<String> oneLine = GFGeneMap.get(i);
//			String gfID_prob= oneLine.get(1);
//
//			for (String s : oneLine) {
//				System.out.print(s + "\t");
//			}
//			
////			for(int k=0; k<oneLine.size()-1;k++){ // to change the old top bottom classifications
////				
////				System.out.print(oneLine.get(k)+"\t");
////			}
////			
////		
////			if(topGFs.contains(gfID_prob)){
////				
////				System.out.print("top"+"\n");
////			}
////			
////			if(bottomGFs.contains(gfID_prob)){
////				
////				System.out.print("bottom"+"\n");
////			}
////			
////			if((!topGFs.contains(gfID_prob))&&(!bottomGFs.contains(gfID_prob))){
////				System.out.print("middle"+"\n");
////			}
//
//			for (int j = 0; j < gfIDs_SC.size(); j++) {
//
//				if (gfIDs_SC.get(j).equals(gfID_prob)) {
//					
//					found =true;
//					
//					System.out.print(groupID.get(j)+"\n");
//
//				}
//			}
//			
//			if(found==false){
//				System.out.print("NaN"+"\n");
//			}
//		}

		
//		String GFcountsFileToSort = "/home/setas/Desktop/setas/Project1/Model_Tests/InputFiles/GFcounts/MGCF5/37spe-MGCF5-9178core-OrderNewickTree";
//		List<List<String>> map = cmmFunct.readMapFile(GFcountsFileToSort);
		
//		for(int i=0; i< gfIDs_SC.size();i++) {
//			
//			for(List<String> lsm: map) {
//				
//				if(lsm.get(0).equalsIgnoreCase(gfIDs_SC.get(i))) {
//					
//					for(int j=0; j<lsm.size();j++) {
//						
//						System.out.print(lsm.get(j)+"\t");
//					}
//					
//					System.out.print("\n");
//				}
//			}
//		}
		
		
		String myFile = "/home/setas/Desktop/GF_GO_All_comma";
		
		List<List<String>> map = cmmFunct.readMapFile(myFile);
		
		for (int i = 0; i < gfIDs_SC.size(); i++) {
			
			System.out.print(gfIDs_SC.get(i));
			
			for (List<String> lsm : map) {

				if (lsm.get(0).equalsIgnoreCase(gfIDs_SC.get(i))) {

					for (int j = 1; j < lsm.size(); j++) {

						System.out.print("\t"+lsm.get(j));
					}
				}
			}
			System.out.print("\n");
		}
		
		
	}
}
