package be.ugent.psb.setas.independent_parsers.RVSEng;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;
import java.util.ArrayList;

public class PrintDecsForMWUGOresults {
	
	public static void main(String[] args) {

	CommonFunctions cmf = new CommonFunctions();

	String path1 = "C:\\Users\\setar\\git\\StochasticBD\\data\\MWU_All_avgPGFD";
	ArrayList<String> GOIDs = cmf.readColX_String(path1, 0);
	ArrayList<String> median_pos = cmf.readColX_String(path1, 1);
	ArrayList<String> median_neg = cmf.readColX_String(path1, 2);
	ArrayList<String> U = cmf.readColX_String(path1, 3);
	ArrayList<String> Z = cmf.readColX_String(path1, 4);
	ArrayList<String> pvalues = cmf.readColX_String(path1, 5);
	ArrayList<String> all = cmf.readColX_String(path1, 6);

	String path2 = "C:\\Users\\setar\\git\\StochasticBD\\data\\ranks_gfID_goHierarch_Desc_inOrderLam";
	ArrayList<String> GOIDs2 = cmf.readColX_String(path2, 2);
    ArrayList<String> Desc2 = cmf.readColX_String(path2, 3);

	for(int i = 0; i< GOIDs.size(); i++){

		boolean found = false;
		String prob = GOIDs.get(i);

		for (int j = 0; j < GOIDs2.size(); j++) {

			if (!found && GOIDs.get(i).equals((GOIDs2).get(j))) {

					System.out.println(prob + "\t" + Desc2.get(j) + "\t" + median_pos.get(i) + "\t" + median_neg.get(i)
							+ "\t" + U.get(i) + "\t" + Z.get(i) + "\t"+pvalues.get(i) + "\t" + all.get(i));
					found = true;
			}
		}
	}
	}
}
