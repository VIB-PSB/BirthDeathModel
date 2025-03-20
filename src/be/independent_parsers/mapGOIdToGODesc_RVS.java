package be.ugent.psb.setas.independent_parsers;

import java.util.ArrayList;

public class mapGOIdToGODesc_RVS {

	public static void main(String[] args) {

		CommonFunctions cmmFunct = new CommonFunctions();

		String path_Go_Desc = "C:\\Users\\setar\\git\\StochasticBD\\data\\ranks_gfID_goHierarch_Desc_inOrderLam";
		ArrayList<String> goIDs_ref = cmmFunct.readColX_String(path_Go_Desc, 2);
		ArrayList<String> dec_ref = cmmFunct.readColX_String(path_Go_Desc, 3);

		String myFile = "C:\\Users\\setar\\git\\StochasticBD\\data\\PGFD_ranks_GO_Bottom";
		ArrayList<String> goIDs = cmmFunct.readColX_String(myFile, 0);
		ArrayList<Double> median_pos = cmmFunct.readColX_double(myFile, 1);
		ArrayList<Double> median_neg = cmmFunct.readColX_double(myFile, 2);
		ArrayList<Double> u = cmmFunct.readColX_double(myFile, 3);
		ArrayList<Double> z = cmmFunct.readColX_double(myFile, 4);
		ArrayList<Double> pvalues = cmmFunct.readColX_double(myFile, 5);
		ArrayList<String> topBottom = cmmFunct.readColX_String(myFile, 6);
		
		ArrayList<String> goIDs_visited = new ArrayList<String>() ;

		for (int i = 0; i < goIDs.size(); i++) {

			String goID_prob = goIDs.get(i);

			for (int j = 0; j < goIDs_ref.size(); j++) {

				if (goID_prob.equals(goIDs_ref.get(j))&& !goIDs_visited.contains(goID_prob)) {
					
					goIDs_visited.add(goID_prob);

					System.out.println(goID_prob + "\t" + dec_ref.get(j) + "\t" + median_pos.get(i) + "\t"
							+ median_neg.get(i) + "\t" + u.get(i) + "\t" + z.get(i) + "\t" + pvalues.get(i) + "\t"
							+ topBottom.get(i));

				}
			}
		}
	}

}
