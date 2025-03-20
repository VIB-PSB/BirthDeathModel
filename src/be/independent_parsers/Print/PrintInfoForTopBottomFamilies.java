package be.ugent.psb.setas.independent_parsers.Print;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class PrintInfoForTopBottomFamilies {

	public static void main(String[] args) {

		CommonFunctions cmmFunct = new CommonFunctions();

		String path1 = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/37AngioSperms/baseFiles/9178coreGF_inOrderLam";
		List<String> GF_id_order = cmmFunct.read1ColFile_String(path1);

		String path = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/TBdupPercent/onlyNonZeroFamilies/GM_nonZero";

		ArrayList<String> gfID = cmmFunct.readColX_String(path, 0);
		ArrayList<Double> isBlock = cmmFunct.readColX_double(path, 1);
		ArrayList<Double> isTandem = cmmFunct.readColX_double(path, 2);
		ArrayList<Double> isBT = cmmFunct.readColX_double(path, 3);

		int i = 8486;
		String gfLastTopGF = GF_id_order.get(i);
//		System.out.println(gfLastTopGF);

//		while (i > 0) { //For Top
		while (i < GF_id_order.size()) {

//			System.out.println(gfLastTopGF);
			
			int index = cmmFunct.searchListString_index(gfLastTopGF, gfID);

			if (index > 0) {

				for (int j = index; j < gfID.size(); j++) {
					System.out.println(gfID.get(j) + "\t" + isBlock.get(j)
							+ "\t" + isTandem.get(j) + "\t" + isBT.get(j));
				
				}
				
				return;
			}

			else {
//				i = i - 1; //For Top families
				i = i+1; // Bottom
				gfLastTopGF = GF_id_order.get(i);

			}

		}

	}

}
