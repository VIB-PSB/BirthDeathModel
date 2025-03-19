package be.ugent.psb.setas.independent_parsers.GO;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class PrintRanksGFGoDescfile {

	/*
	 * read first line: print (i=)1 \t line read second line, if GF-ID =
	 * GF-ID-previous, type i \t line else print i+1 \t line
	 */

	public double findRankOnOtherFile(String rankFile, String GFid_probe){
		
		CommonFunctions cmf = new CommonFunctions();
		
		double foundNewRank = 0;
		
		ArrayList<String> gfIDs_otherFile = cmf.readColX_String(rankFile, 0);
		ArrayList<Double> newRanks = cmf.readColX_double(rankFile, 1);
		
		for(int k=0; k<gfIDs_otherFile.size();k++){
			
			if(gfIDs_otherFile.get(k).equals(GFid_probe)){
			
				foundNewRank = newRanks.get(k);
			}
			
		}
		
		return foundNewRank;
		
	}
	
	
	public static void main(String[] args) {

		CommonFunctions cmmFunc = new CommonFunctions();
		PrintRanksGFGoDescfile prntrank = new PrintRanksGFGoDescfile();

		// you assume that this file is in correct order according to lambdas +
		// all columns have the same length

		// String path
		// ="/home/setas/git/IndependentParsers/Independent Parsers/src/files/GF_GO_Desc_allHierarchy_28EudGFs_inOrder";
//		String path = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/GFidOrderLam_GOid_Desc_GoSlimPlants_37spe";
//		String path="/home/setas/git/IndependentParsers/Independent Parsers/src/files/GFidOrderLam_GOid_Desc_GoSlimPlants_28Eud";
		
		String path="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/MannWhitneyUtest/new1000TopBottom_combRankBlockLambda/ranksOld_gfID_goHierarch_Desc_inOrderLam";
		
		ArrayList<String> gfIDs = cmmFunc.readColX_String(path, 1);
		ArrayList<Integer> goIDs = cmmFunc.readColX_Int(path, 2);
		ArrayList<String> goDesc = cmmFunc.readColX_String(path, 3);

		// System.out.println("size"+gfIDs.size());

		int rank = 1;
		System.out.println(rank + "\t" + gfIDs.get(0) + "\t" + goIDs.get(0)
				+ "\t" + goDesc.get(0));

//		for (int i = 1; i < gfIDs.size(); i++) {
//
//			if (gfIDs.get(i).equals(gfIDs.get(i - 1))) {
//
//				System.out.println(rank + "\t" + gfIDs.get(i) + "\t"
//						+ goIDs.get(i) + "\t" + goDesc.get(i));
//			}
//
//			else {
//				rank = rank + 1;
//				System.out.println(rank + "\t" + gfIDs.get(i) + "\t"
//						+ goIDs.get(i) + "\t" + goDesc.get(i));
//			}
//		}
		
		String rankFile="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/MannWhitneyUtest/new1000TopBottom_combRankBlockLambda/GFid_avgRankLambdaBlock";
		
		for (int i = 1; i < gfIDs.size(); i++) {
			
			double newRank = prntrank.findRankOnOtherFile(rankFile,gfIDs.get(i));
			
			System.out.println(newRank+"\t"+gfIDs.get(i)+"\t"+goIDs.get(i)+"\t"+goDesc.get(i));
			
		}

	}
}