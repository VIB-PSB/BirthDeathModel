package be.ugent.psb.setas.independent_parsers;

import java.util.ArrayList;
import java.util.List;

public class CalcAvgDiffWithIdealProfiles {
	
	
//	public static List<Integer> absArrayList(List<Integer> myArrayList){
//		
////		System.out.println(myArrayList.size());
//		
//		List<Integer> absarr = new ArrayList<Integer>(myArrayList.size());
//		
////		System.out.println("size"+absarr.size()); //size =0 Problem
//		
//		for(int i=0; i<myArrayList.size();i++){
//			
//			int absvalue = Math.abs(myArrayList.get(i));
//			
//			absarr.add(absvalue);
//		}	
//		return absarr;		
//	}

	public static int[] convertArryListToArray(ArrayList<Integer> a){
	
		int [] arrayForm = new int[a.size()];
		
		for(int i=0; i<a.size(); i++){	
			arrayForm[i] = a.get(i).intValue();
		}
		return arrayForm;
	}
	
	public static int[] convertArryListToArray_2(List<Integer> a){	
		
		int [] arrayForm = new int[a.size()];
		
		for(int i=0; i<a.size(); i++){
			arrayForm[i] = a.get(i).intValue();
		}
		return arrayForm;
	}

//	public static void main(String[] args) {
//
//		CommonFunctions cmmFunct = new CommonFunctions();
//		GFsizeDistribution gfsizedist = new GFsizeDistribution();
//
//		String path1 = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/37AngioSperms/baseFiles/9178coreGF_inOrderLam";
////		String path1 = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/GFsOrder_37spe_Tau_Mon2_Musa3";
//		List<String> GF_ids_order = cmmFunct.read1ColFile_String(path1);
//
////		String path2 = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/combinedOutput/comSortLam_37speMGCF5_CPMpval_9178coreGF";
//		String path2= "/home/setas/Desktop/setas/Project1/Results/CompareRankings/newCombinedOutput/combinedOutputFiles_rawSorted/SortedLambda/37spe/comOut_37spe_Tau_Mon2_Musa3_SortedLambda";
//
//		ArrayList<String> gfIDs_combinedOutput = cmmFunct.readColX_String(
//				path2, 0);
//		ArrayList<Integer> optimalRootSizes = cmmFunct.readColX_Int(path2, 1);
//
//		String idealProfilePath = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ComparisonsWithIdealprofiles/IdealProfiles_R1";
////
////		ArrayList<String> speInOrder = cmmFunct.readColX_String(
////				idealProfilePath, 0);
////		ArrayList<Integer> idealProfileR1 = cmmFunct.readColX_Int(
////				idealProfilePath, 1);
//		ArrayList<Integer> idealProfileR1 = cmmFunct.readColX_Int(
//				idealProfilePath, 4); // Tau-Mon2_Musa3
//
//		ReadGFcountsFile rgf = new ReadGFcountsFile();
//		List<List<Integer>> gfCounts = rgf
//				.read_all("/home/setas/git/StochasticBD/src/Files/GeneFamilyCounts/37spe-MGCF5-9178core-OrderNewickTree");
//		ArrayList<String> gfIDs_geneCountFile = rgf.getGfIDs();
//		
//		DecimalFormat df = new DecimalFormat("0.0000");
//
//		for (int i = 0; i < GF_ids_order.size(); i++) { // to be able to have
//														// all output files in
//														// one order
//			String gfID_current = GF_ids_order.get(i);
//
//			int indexInComOutFile = cmmFunct.searchListString_index(
//					gfID_current, gfIDs_combinedOutput);
//			int rootSize_current = optimalRootSizes.get(indexInComOutFile);		
//			System.out.print(gfID_current+"\t");
//			
//			ArrayList<Integer> idealProfile_multiplyByRstar = gfsizedist.multiplyByXList(idealProfileR1,
//					rootSize_current);
//			
//			int [] idealProf = convertArryListToArray(idealProfile_multiplyByRstar);
//
//			for (int j = 0; j < gfIDs_geneCountFile.size(); j++) {
////				if (gfIDs_geneCountFile.get(j).equalsIgnoreCase(gfID_current)) {
////					List<Integer> geneCountProfile_current = gfCounts.get(j);
////
////					ArrayList<Integer> diffIdealOptR = gfsizedist
////							.getDiffTwoIntList((ArrayList<Integer>) geneCountProfile_current,idealProfile_multiplyByRstar);
////					
////					int [] diff = convertArryListToArray(diffIdealOptR);
////					int [] absdiff = new int[diff.length];
////					int sumAbs_stmae =0;
////					for(int k=0; k<37;k++){
////						absdiff[k]= Math.abs(diff[k]);
////					}
////					double [] normalized = new double[absdiff.length];
////					double sum=0;
////					for(int l=0; l<37;l++){	
////						normalized[l]= absdiff[l]*1.000/idealProf[l]*1.000;	
////						sum+=normalized[l];
//////						System.out.print(df.format(normalized[l])+"\t");
////					}
////					System.out.print(df.format(sum*1.0000/37.0000));
////					System.out.println();
////				}
//			}
//		}
//
//	}

}
