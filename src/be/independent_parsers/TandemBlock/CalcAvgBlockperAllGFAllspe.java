package be.ugent.psb.setas.independent_parsers.TandemBlock;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class CalcAvgBlockperAllGFAllspe {

	public static void main(String[] args) {

		CommonFunctions cmmFunct = new CommonFunctions();

		String path1 = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/37AngioSperms/baseFiles/9178coreGF_inOrderLam";
		List<String> GF_id_order = cmmFunct.read1ColFile_String(path1);

		 String path =
		 "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/new/isTBAllGenesAllGFs_37spePlazaID";
//		String path = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/NewFits/LoessOutput2.txt";

		ArrayList<String> gfID = cmmFunct.readColX_String(path, 0);
		ArrayList<String> geneID = cmmFunct.readColX_String(path, 1);
		ArrayList<Integer> isBlock = cmmFunct.readColX_Int(path, 2);
		ArrayList<Integer> isTandem = cmmFunct.readColX_Int(path, 3);
		ArrayList<Integer> isTB = cmmFunct.readColX_Int(path, 4);
		
//		ArrayList<Double> Residue = cmmFunct.readColX_double(path, 6); // in case of exp div. it is equal to the number of pairs

		DecimalFormat df = new DecimalFormat("0.000");

		for (int j = 0; j < GF_id_order.size(); j++) {

			String gfcurrent = GF_id_order.get(j);
//			int totalNumDuplicates = 0;
			int totalNumGenesInBlockDup = 0;
			int totalNumGenesInTandemDup = 0;
			int totalNumGenesInBTDup = 0;

			int totalNumGenesInGF = 0;
			
//			double sumAbsExpDiv =0;

			for (int i = 0; i < gfID.size(); i++) {

				if (gfID.get(i).equals(gfcurrent)) { // The same family

					totalNumGenesInGF += 1;
					
//					sumAbsExpDiv+= Residue.get(i);
					

//					// if (isBlock.get(i) == 1.0 || isTandem.get(i) == 1 ||
//					// isTB.get(i)==1) {
//					// totalNumDuplicates = totalNumDuplicates + 1; // we only
//					// count the ones that are some sort of duplicate in the
//					// calculation
//					//
//					// if(isBlock.get(i)==1){
//					// totalNumGenesInBlockDup +=1;
//					// }
//					// }
//
					if (isBlock.get(i) == 1 && isTandem.get(i) == 0
							&& isTB.get(i) == 0) {

						totalNumGenesInBlockDup += 1;

					}

					if (isBlock.get(i) == 0 && isTandem.get(i) == 1
							&& isTB.get(i) == 0) {

						totalNumGenesInTandemDup += 1;

					}

					if (isBlock.get(i) == 1 && isTandem.get(i) == 1
							&& isTB.get(i) == 1) {

						totalNumGenesInBTDup += 1;

					}

				}

			}

			// double avgBlockDup = totalNumGenesInBlockDup * 1.000/
			// totalNumDuplicates * 1.000;
			// System.out.println(gfcurrent+"\t"+totalNumGenesInBlockDup+"\t"+totalNumDuplicates);
			// System.out.println(gfcurrent + "\t" + df.format(avgBlockDup));

			double avgBlockDup = totalNumGenesInBlockDup * 1.000
					/ totalNumGenesInGF * 1.000;
			double avgTandemDup = totalNumGenesInTandemDup * 1.000
					/ totalNumGenesInGF * 1.000;
			double avgBTDup = totalNumGenesInBTDup * 1.000
					/ totalNumGenesInGF * 1.000;

			System.out.println(gfcurrent + "\t" + df.format(avgBlockDup) + "\t"
					+ df.format(avgTandemDup) + "\t" + df.format(avgBTDup));
			
			
//			if(totalNumGenesInGF!=0){
			
//			double avgExpDiv = sumAbsExpDiv*1.000 / totalNumGenesInGF*1.000;

//			System.out.println(gfcurrent+"\t"+totalNumGenesInGF+"\t"+df.format(avgExpDiv));
//			}
			
//			else{
//				System.out.println(gfcurrent+"\t"+totalNumGenesInGF+"\t"+"N/A");
//			}
		}

	}

}
