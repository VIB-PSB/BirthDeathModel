package be.ugent.psb.setas.independent_parsers.TandemBlock;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class TBpercentagesAllGFs {

	public static void main(String[] args) {

		CommonFunctions cmmFunct = new CommonFunctions();

		String path1 = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/37AngioSperms/baseFiles/9178coreGF_inOrderLam";
		List<String> GF_id_order = cmmFunct.read1ColFile_String(path1);

//		String path = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/TBpercent_12Spe_AllGFsInOrder";
		
//		String path = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/mibel-newIAdhore/ProcessedFiles/TB_Top709_29spePlaza";
//		String path="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/mibel-newIAdhore/ProcessedFiles/TB_Bottom692_29spePlaza";
//		String path="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/Tpar/TndemBlock_29spePlaza_Tpar_Top709";
//		String path="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/new/AllGenes_AllGFs_37pePlazaIDs";
		String path="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/new/isTB_7RemainedSpe_AllGFs";

		ArrayList<String> gfID = cmmFunct.readColX_String(path,0);
		ArrayList<String> geneID = cmmFunct.readColX_String(path,1);
		ArrayList<Double> isBlock = cmmFunct.readColX_double(path,2);
		ArrayList<Double> isTandem = cmmFunct.readColX_double(path,3);
		ArrayList<Double> isTB = cmmFunct.readColX_double(path,4);
		
		String geneIDStarter = "CM";

		for (int j = 0; j < GF_id_order.size(); j++) {
//		for (int j = 0; j < 709; j++) { //for the top families
//		for (int j = 8684; j < GF_id_order.size(); j++) { //for the bottom families
		
			String gfcurrent = GF_id_order.get(j);
			
//			if(gfcurrent.split("_").length > 1){ //this way you are ignoring the differences of all _1 and _2 gene families				
//				String [] gfsOrder = GF_id_order.get(j).split("_");
//				gfcurrent = gfsOrder[0];	}

			int isB = 0;
			int isT = 0;
			int isBT = 0;
			int totalNumDuplicates=0;
			int totalNumGenesInSpe =0;

			for (int i = 0; i < gfID.size(); i++) {
				
				if (gfID.get(i).equals(gfcurrent)) { // The same family
					
//					if(geneID.get(i).split("ATR").length <= 1){ // only for Athaliana to remove the Amborella genes: ATR

					if (geneID.get(i).split(geneIDStarter).length > 1) { 
						
						totalNumGenesInSpe = totalNumGenesInSpe + 1;  // we changed it later to calculate the percentages only based on known duplicates

						if (isBlock.get(i) == 1.0 && isTandem.get(i) == 1) {

							totalNumDuplicates = totalNumDuplicates + 1; // we only count the ones that are some sort of duplicate in the calculation
							isBT = isBT + 1;

						}

						if (isBlock.get(i) == 1 && isTandem.get(i) == 0) {

							totalNumDuplicates = totalNumDuplicates + 1;
							isB = isB + 1;
						}

						else if (isBlock.get(i) == 0 && isTandem.get(i) == 1) {

							totalNumDuplicates = totalNumDuplicates + 1;
							isT = isT + 1;
						}

//					}
					}
				}
				
				
			}

			if(totalNumDuplicates != 0){ //to print the large combined file, we need to include 0,0,0 to have it for all families
				
			DecimalFormat df = new DecimalFormat("0.000");		
//			double isBpercent = isB*1.00/totalNumDuplicates;
//			double isTpercent = isT*1.00/totalNumDuplicates;
//			double isBTpercent = isBT*1.00/totalNumDuplicates;			
//			System.out.println(gfcurrent + "\t" + df.format(isBpercent) + "\t" + df.format(isTpercent)+"\t"+ df.format(isBTpercent));
			
			
			// for a combined score for all spe later:	
			System.out.println(gfcurrent + "\t" + isB + "\t" + totalNumDuplicates);
			
			}
			else if(totalNumDuplicates == 0 && totalNumGenesInSpe !=0){
//				System.out.println(gfcurrent + "\t" + isB + "\t" + isT + "\t"
//						+ isBT);
				
				System.out.println(gfcurrent + "\t" + isB+"\t"+ totalNumDuplicates);
			}
			
			else if(totalNumDuplicates == 0 && totalNumGenesInSpe ==0){
//				System.out.println(gfcurrent + "\t" + "NA" + "\t" + "NA" + "\t"
//						+ "NA");
				
				System.out.println(gfcurrent + "\t" + "NA"+"\t"+"NA");
			}
						
		}

	}

}
