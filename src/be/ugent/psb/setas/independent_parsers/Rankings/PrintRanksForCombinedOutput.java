package be.ugent.psb.setas.independent_parsers.Rankings;

import java.util.ArrayList;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class PrintRanksForCombinedOutput {

	public static void main(String[] args) {

		CommonFunctions cmmFunc = new CommonFunctions();
		
		/** The combined output files must be sorted in order of lambda **/
//		String path="/home/setas/git/IndependentParsers/Independent Parsers/src/files/comOut_28Eud_Prun37_9178OrthGF_SortLam";
		String path="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-8Monocots/9178_8Monocots_sorted_lambda";
//		String path="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-28Eudicots/combinedOutput/comOut_28Eud_Prun37_9178OrthGF_SortLam";
//		String path="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/combinedOutput/comSortLam_37speMGCF5_CPMpval_9178coreGF";
		
		ArrayList<String> gfIDs = cmmFunc.readColX_String(path,0);
		ArrayList<Integer> rootSize = cmmFunc.readColX_Int(path,1);
		ArrayList<Double> lambdas = cmmFunc.readColX_double(path,2);
		ArrayList<Double> likelihoods = cmmFunc.readColX_double(path,3);
		ArrayList<Double> pvalues = cmmFunc.readColX_double(path,4);

		int rank = 1;
		int goDownInRank =0;
		
		System.out.println(rank + "\t" + gfIDs.get(0) + "\t" + rootSize.get(0)
				+ "\t" + lambdas.get(0) + "\t" + likelihoods.get(0) + "\t"
				+ pvalues.get(0));

		for (int i = 1; i < gfIDs.size(); i++) { // the same gene family, pritn the exact same rank

			if (gfIDs.get(i).equals(gfIDs.get(i - 1)) ) {

				System.out.println(rank + "\t" + gfIDs.get(i) + "\t"
						+ rootSize.get(i) + "\t" + lambdas.get(i) + "\t"
						+ likelihoods.get(i) + "\t" + pvalues.get(i));
			}
			

			else if(lambdas.get(i)-lambdas.get(i-1) < 0.000001){ //if lambda is the same, calculate an average rank
				
//				System.out.println("lambdas equal");
				goDownInRank +=1;	
				System.out.println(rank + "\t" + gfIDs.get(i) + "\t"
						+ rootSize.get(i) + "\t" + lambdas.get(i) + "\t"
						+ likelihoods.get(i) + "\t" + pvalues.get(i));
			}
			else {

				goDownInRank +=1;
				rank = 1 + goDownInRank;
				System.out.println(rank + "\t" + gfIDs.get(i) + "\t"
						+ rootSize.get(i) + "\t" + lambdas.get(i) + "\t"
						+ likelihoods.get(i) + "\t" + pvalues.get(i));
			}
		}
	}
}
