package be.ugent.psb.setas.independent_parsers;

import java.util.ArrayList;
import java.util.List;

public class MapAthaGenes {
	
	public static boolean searchArrayOfString(String s, String[] myArray){
		
		boolean f =false;
		
		for(int j=0; j< myArray.length;j++){
			
			if(s.equalsIgnoreCase(myArray[j])){
				
				f=true;
			}
		}
		
		return f;
	}

	public static void main(String[] args) {
		CommonFunctions cmmFunct = new CommonFunctions();

		// String path_GF_Ath =
		// "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/GFid_AthaGene_InOrder";
		//
		// String path_TF_Ath_family =
		// "/home/setas/Desktop/setas/setas_stmae_locar/Reviews/Ath_TF_list.txt";
		//
		// String path_TF_family_retainedFreq =
		// "/home/setas/Desktop/setas/setas_stmae_locar/Reviews/TF_family_retainedFreq_freeling";
		//
		// List<List<String>> map_GF_Ath = cmmFunct.readMapFile(path_GF_Ath);
		// List<List<String>> map_TF_Ath_family =
		// cmmFunct.readMapFile(path_TF_Ath_family);
		// List<List<String>> map_TF_family_retainedFreq =
		// cmmFunct.readMapFile(path_TF_family_retainedFreq);
		//
		// for (int i = 0; i < map_GF_Ath.size(); i++) {
		//
		// List<String> GF_Ath = map_GF_Ath.get(i);
		//
		// String GFid = GF_Ath.get(0);
		// String Atha = GF_Ath.get(1);
		// boolean found = false;
		// boolean hasRetentionRate = false;
		//
		// for (List<String> ls_TF_Ath_family : map_TF_Ath_family) {
		//
		// if (ls_TF_Ath_family.get(1).equalsIgnoreCase(Atha)) {
		//
		// found = true;
		// System.out.print(GFid + "\t" + Atha + "\t" +
		// ls_TF_Ath_family.get(2));
		//
		// for (List<String> ls_TF_retainedFreq : map_TF_family_retainedFreq) {
		//
		// if
		// (ls_TF_retainedFreq.get(0).equalsIgnoreCase(ls_TF_Ath_family.get(2)))
		// {
		//
		// hasRetentionRate = true;
		// System.out.print("\t" + ls_TF_retainedFreq.get(1));
		// }
		//
		// }
		//
		// if (hasRetentionRate == false) {
		// System.out.print("\t" + "NA");
		// }
		//
		// System.out.print("\n");
		// }
		// }
		//
		// if (found == false) {
		// System.out.println(GFid + "\t" + Atha + "\t" + "NA" + "\t" + "NA");
		// }
		// }

		// String path_GFinOrder =
		// "/home/setas/Desktop/setas/setas_stmae_locar/Reviews/compareLi2015/AllGFs_InOrderLam_MusaFirst2Close_PartitionSC";
		// String path_GF_Ath_TFfamily_retention =
		// "/home/setas/Desktop/setas/setas_stmae_locar/Reviews/Freeling/GF_Ath_TF_family_retentionFreq";
		//
		// List<List<String>> map_GF_order =
		// cmmFunct.readMapFile(path_GFinOrder);
		// List<List<String>> map_TF_Ath_TFfamily_retention =
		// cmmFunct.readMapFile(path_GF_Ath_TFfamily_retention);
		//
		//
		// for (int i = 0; i < map_GF_order.size(); i++) {
		//
		// String gfID = map_GF_order.get(i).get(0);
		// boolean found =false;
		//
		// for (List<String> ls : map_TF_Ath_TFfamily_retention) {
		//
		// if (ls.get(0).equalsIgnoreCase(gfID)) {
		//
		// found = true;
		//
		// for (String s : ls) {
		// System.out.print(s + "\t");
		// }
		// System.out.print("\n");
		// }
		// }
		// if(found==false){
		//
		// System.out.println(gfID+"\t"+"NA"+"\t"+"NA"+"\t"+"NA");
		// }
		// }

		// MapAthaGenes mapath = new MapAthaGenes();
		// String path_GF_Ath_TFF_retention_order =
		// "/home/setas/Desktop/setas/setas_stmae_locar/Reviews/Freeling/ALLGF_Ath_TF_family_retentionFreq_inOrder";
		// List<List<String>> map_GF_order =
		// cmmFunct.readMapFile(path_GF_Ath_TFF_retention_order);
		//
		// for (int counter = 0; counter < 9178; counter++) {
		//
		// List<String> ls = map_GF_order.get(counter);
		//
		// String gfID = ls.get(0);
		// String Ath = ls.get(1);
		// String TFfamily = ls.get(2);
		// String retention = null;
		//
		// if (ls.get(3) != null) {
		// retention = ls.get(3);
		// }
		//
		// System.out.print(gfID + "\t" + Ath);
		//
		// while (map_GF_order.get(counter + 1).get(0).equalsIgnoreCase(gfID)) {
		//
		// System.out.print("," + map_GF_order.get(counter + 1).get(1));
		// counter = counter + 1;
		//
		// }
		//
		// System.out.print("\t" + TFfamily + "\t" + retention + "\n");
		// }

		String path_com_Output = "/home/setas/Desktop/setas/paper-R5/MBE_submission/Supplementary_Tables/tableS1_combOutput_MusaFirst2_noHeader";
//		String path_TF_Ath_family = "/home/setas/Desktop/setas/setas_stmae_locar/Reviews/Ath_TF_list.txt";
//		String path_TF_family_retainedFreq = "/home/setas/Desktop/setas/setas_stmae_locar/Reviews/TF_family_retainedFreq_freeling";

		List<List<String>> map = cmmFunct.readMapFile(path_com_Output);
//		List<List<String>> map_TF_Ath_family = cmmFunct.readMapFile(path_TF_Ath_family);
//		List<List<String>> map_TF_family_retainedFreq = cmmFunct.readMapFile(path_TF_family_retainedFreq);

//		for (List<String> ls : map) {
//
//			String[] ath_genes = ls.get(3).split(",");
//			ArrayList<String> TFFs = new ArrayList<String>();
//		
//
//			for (int j = 0; j <= 3; j++) {
//				System.out.print(ls.get(j) + "\t");
//			}
//			
//			System.out.print(ls.get(7)+"\t"+ls.get(9)+"\t");
//
//			for (String s : ath_genes) {
//
//				for (List<String> TF_family : map_TF_Ath_family) {
//
//					if (TF_family.get(1).equalsIgnoreCase(s) && (!TFFs.contains(TF_family.get(2)))) {
//
//						TFFs.add(TF_family.get(2));
//					}
//				}
//			}
//			
//			if(TFFs.size()==0){
//				TFFs.add("NA");
//			}
//			for (String ars : TFFs) {
//				
//				boolean retentionFound=false;
//
//				System.out.print(ars + "\t");
//
//				for (List<String> retentions : map_TF_family_retainedFreq) {
//
//					if (retentions.get(0).equals(ars)) {
//						
//						retentionFound=true;
//						System.out.print(retentions.get(1));
//					}
//				}
//				
//				if(retentionFound==false){
//					System.out.print("NA");
//				}
//			}
//
//			System.out.print("\n");
//		}
		
	
		
		String path_TFfamily="/home/setas/Desktop/setas/setas_stmae/Reviews/Freeling/Freeling2008/Freeling 2008";
		List<List<String>> map_TFfamily = cmmFunct.readMapFile(path_TFfamily);
		
        
		for (int i = 0; i < map.size(); i++) {

			List<String> combOutputLine = map.get(i);
			String athaGenes = combOutputLine.get(3);
			String[] athaGenes_array = athaGenes.split(",");
			
			List<String> families = new ArrayList<String>();
			boolean noFamiliesFound = true;

			for (List<String> lss : map_TFfamily) {

				
				if(lss.size()>=12){
					
					if(lss.get(11).trim().length() >0){
					
//					System.out.println(combOutputLine.get(0)+"\t"+lss.get(11));
					if (searchArrayOfString(lss.get(1),athaGenes_array) && (!families.contains(lss.get(11)))) { //1st item is the atha genes

						  noFamiliesFound = false;
	                      families.add(lss.get(11));
					}
					}
				}
				
			}

			for(String s:combOutputLine){
			System.out.print(s+"\t");
			}
			
			for(String f:families){
				System.out.print(f+"\t");
			}
			
			if(noFamiliesFound==true){
				System.out.print("NA"+"\t");
			}
			
			System.out.print("\n");
		}

	}

}
