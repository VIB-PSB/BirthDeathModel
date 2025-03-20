package be.ugent.psb.setas.independent_parsers.RVSEng;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class RankTheRvsEngAnalysedOutput {

	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();

		String path1 = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/37AngioSperms/baseFiles/9178coreGF_inOrderLam";

		ArrayList<String> gfIDsOrder = cmf.read1ColFile_String(path1);

		String path2 = "/home/setas/Desktop/setas/Project1/Results/RvsEng/Removing1WGD/LamBiggerLoglkLower/LambdaBiggerLogLkLower_rmWGD17";

		ArrayList<String> gfIDs = cmf.readColX_String(path2, 0);
		ArrayList<Double> rStar = cmf.readColX_double(path2, 1);
		ArrayList<Double> lambdaStar = cmf.readColX_double(path2, 2);
		ArrayList<Double> loglkStar = cmf.readColX_double(path2, 3);
		ArrayList<Double> loglkRemoval = cmf.readColX_double(path2, 4);
		ArrayList<String> comparison = cmf.readColX_String(path2, 5);
		ArrayList<Double> rStarNew = cmf.readColX_double(path2, 6);
		ArrayList<Double> lambdaStarNew = cmf.readColX_double(path2, 7);
		ArrayList<Double> loglkStarNew = cmf.readColX_double(path2, 8);
		ArrayList<String> comparisonNew = cmf.readColX_String(path2, 3);

		for (int i = 0; i < gfIDsOrder.size(); i++) { //i = real rank

			String gfId_prob = gfIDsOrder.get(i);

			if (cmf.searchListString_boolean(gfId_prob, gfIDs)) {

				int index = cmf.searchListString_index(gfId_prob, gfIDs);

				System.out.println(i + "\t" + gfId_prob + "\t"
						+ rStar.get(index) + "\t" + lambdaStar.get(index)
						+ "\t" + loglkStar.get(index) + "\t"
						+ loglkRemoval.get(index) + "\t"
						+ comparison.get(index) + "\t" + rStarNew.get(index)
						+ "\t" + lambdaStarNew.get(index) + "\t"
						+ loglkStarNew.get(index) + "\t"
						+ comparisonNew.get(index));

			}

		}

	}
}
