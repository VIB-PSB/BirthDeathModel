package be.ugent.psb.setas.independent_parsers.Rankings;

import java.util.ArrayList;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class CompareRankAngio9178_28Eudicot {

	public static void main(String[] args) {

		CommonFunctions cmmFunc = new CommonFunctions();
		ArrayList<String> gfIDs_original = cmmFunc
				.read1ColFile_String("/home/setas/git/IndependentParsers/Independent Parsers/src/files/9178coreGF_inOrderLam");

		String path = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/ranks_comOut_28Eud_Prun37_9178OrthGF_SortLam";
		ArrayList<Integer> ranks = cmmFunc.readColX_Int(path, 0);
		ArrayList<String> gfIDs = cmmFunc.readColX_String(path, 1);
		ArrayList<Integer> rootSize = cmmFunc.readColX_Int(path, 2);
		ArrayList<Double> lambdas = cmmFunc.readColX_double(path, 3);
		ArrayList<Double> likelihoods = cmmFunc.readColX_double(path, 4);
		ArrayList<Double> pvalues = cmmFunc.readColX_double(path, 5);

		for (int i = 0; i < gfIDs_original.size(); i++) {

			String gf = gfIDs_original.get(i);

			for (int j = 0; j < gfIDs.size(); j++) {

				if (gfIDs.get(j).equals(gf)) {

					System.out.println((i + 1) + "\t" + ranks.get(j));

				}
			}

		}
	}

}
