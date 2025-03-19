package be.ugent.psb.setas.independent_parsers.Sort;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class SortBasedOnGFsOrderLambda {

	public static void main(String[] args) {

		CommonFunctions cmmFunc = new CommonFunctions();

//		String path1 = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/37AngioSperms/baseFiles/9178coreGF_inOrderLam";
		String path1 ="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/GFsOrder_37spe_Tau_Mon2_Musa3";
//		String path1 = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/GFsOrder_28Eud_37spe_9178coreGF";
		List<String> GFid_order = cmmFunc.read1ColFile_String(path1);
		
//		System.out.println(GFid_order.size());

//		String path2 = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/GFid_GOid_Desc_GoSlimPlants_37spe";
//		String path2="/home/setas/git/IndependentParsers/Independent Parsers/src/files/28EudGoSlimPlantObo_GF_GO_Desc";
//		ArrayList<String> gfIDs = cmmFunc.readColX_String(path2, 0);
//		ArrayList<Integer> goIDs = cmmFunc.readColX_Int(path2, 1);
//		ArrayList<String> goDesc = cmmFunc.readColX_String(path2, 2);
		
		
//		String path2="/home/setas/Desktop/printRanksAvgDistToIdeal";
//		ArrayList<Integer> ranks = cmmFunc.readColX_Int(path2, 0);
//		ArrayList<String> gfIDs = cmmFunc.readColX_String(path2, 1);
//		ArrayList<Double> blockPrct = cmmFunc.readColX_double(path2, 2);
		
//		for (int j = 0; j < 9178; j++) {
//
//			String gfId_prob = GFid_order.get(j);
//
////			if (cmmFunc.searchListString_boolean(gfId_prob, gfIDs)) {
//
//				for (int i = 0; i < gfIDs.size(); i++) {
//
//					if (gfIDs.get(i).equals(gfId_prob)) {
//
//						System.out.print(gfId_prob + "\t" + ranks.get(i) + "\t"
//								+ blockPrct.get(i));
//						System.out.println();
//					}
//				}
//
////			}
//
//		}
		

		
//		String path2="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier2-8Eudicots/ComOut_8EudicotsOnly";
//		ArrayList<String> gfIDs = cmmFunc.readColX_String(path2, 0);
//		ArrayList<Integer> rootSize = cmmFunc.readColX_Int(path2, 1);
//		ArrayList<Double> lambdas = cmmFunc.readColX_double(path2, 2);
//		ArrayList<Double> lks = cmmFunc.readColX_double(path2, 3);
//		ArrayList<Double> pvs = cmmFunc.readColX_double(path2, 4);
//		
//		for (int j = 0; j < 9178; j++) {
//
//			String gfId_prob = GFid_order.get(j);
//
//				for (int i = 0; i < gfIDs.size(); i++) {
//
//					if (gfIDs.get(i).equals(gfId_prob)) {
//
//						System.out.print(gfId_prob + "\t" + rootSize.get(i) + "\t"
//								+ lambdas.get(i)+"\t"+lks.get(i)+"\t"+pvs.get(i));
////						System.out.print(gfId_prob + "\t" + rootSize.get(i) + "\t"
////								+ lambdas.get(i)+"\t"+lks.get(i));
//						System.out.println();
//					}
//				}
//
//		}
		
//		String path2="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/KnKs/referenceLine/allGFs_all37spe_KnKsResidues.txt";
//		ArrayList<String> gfIDs = cmmFunc.readColX_String(path2, 0);
//		ArrayList<Double> avgResidues = cmmFunc.readColX_double(path2, 1);
//
//		for (int j = 0; j < 9178; j++) {
//
//			String gfId_prob = GFid_order.get(j);
//			boolean found = false;
//
//				for (int i = 0; i < gfIDs.size(); i++) {
//					if (gfIDs.get(i).equals(gfId_prob)) {
//						found = true;
//						System.out.println(gfId_prob + "\t" + avgResidues.get(i));
//					}	
//				}
//				
//				if (! found){
//					System.out.println(gfId_prob + "\t" + "NA");
//				}
//
//		}
		
//		String path2= args[0];
//		String path2="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/KnKs/GFsOrder_AllspeGenePairs/Ks5_minKs/AlyrPairs_Ks5";
		String path2 ="/home/setas/Desktop/setas/Project1/RvsEng/Musa3TauMon2/OutputMerged/output_merged_WGD16";	
		List<List<String>> map = cmmFunc.readMapFile(path2);
		ArrayList<String> gfIDs = cmmFunc.readColX_String(path2, 0);
		
//		System.out.println(gfIDs.size());
		
		
		for(int j=0; j<9178; j++){
//		for (int j = 0; j < 1122; j++) { // For the top ones only new default
//		for (int j = 8209; j < 9178; j++) { // For the top ones only new default

			String gfId_prob = GFid_order.get(j);
//			boolean found = false;
			
			for (int i = 0; i < gfIDs.size(); i++) {
				
				if (gfIDs.get(i).equals(gfId_prob)) {
					
					System.out.print(gfId_prob);
					
					for(int k=1; k<map.get(i).size();k++){ // becasue the GF id is printed before start from 1
						
					System.out.print("\t"+map.get(i).get(k));
					
					}
					System.out.print("\n");
					
//					found = true;
				}
			}
			
//			if(found==false){
//				System.out.print(gfId_prob+"\n");	
//			}
			
		}
	}
}
