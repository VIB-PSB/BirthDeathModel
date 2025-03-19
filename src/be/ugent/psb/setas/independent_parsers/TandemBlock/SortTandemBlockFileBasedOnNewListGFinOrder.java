package be.ugent.psb.setas.independent_parsers.TandemBlock;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class SortTandemBlockFileBasedOnNewListGFinOrder {

	public static void main(String[] args) {
		if (args.length != 3){
			System.err.println("Usage: sort_tandem_block GENE_FAMILY_LIST SPECIES_TANDEM_BLOCK NUMBER_OF_TOP");
			return;
		}
		
		CommonFunctions cmmFunc = new CommonFunctions();

//		String path1 = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/37AngioSperms/baseFiles/9178coreGF_inOrderLam";
		String path1 = args[0];
		List<String> GFid_order = cmmFunc.read1ColFile_String(path1);

//		String path = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/new/AllGFs_including0Counts/AllTBpercent_Inc0/AL-allGFs-Include0_TBpercent";
		String path = args[1];
		
		ArrayList<String> gfIDs = cmmFunc.readColX_String(path, 0);
		ArrayList<String> blockP = cmmFunc.readColX_String(path, 1);
		ArrayList<String> tandemP = cmmFunc.readColX_String(path, 2);
		ArrayList<String> tbP = cmmFunc.readColX_String(path, 3);
		
//		int lastTopFamily = Integer.parseInt(args[2]);

		int firstBottomFamily = Integer.parseInt(args[2]);
		
//		for (int j = 0; j < lastTopFamily; j++) { // for the top families
		for (int j = firstBottomFamily; j < 9178; j++) { // for the bottom families

			String gfId_prob = GFid_order.get(j);

			for (int i = 0; i < gfIDs.size(); i++) {

				if (gfIDs.get(i).equals(gfId_prob)) {
					
					if ( !(blockP.get(i).equals("0") ||  blockP.get(i).equals("NA"))) {
						System.out.print(gfId_prob + "\t" + blockP.get(i) + "\t"
								+ tandemP.get(i) + "\t" + tbP.get(i));
						System.out.println();
					}
				}
			}

		}

	}

}
