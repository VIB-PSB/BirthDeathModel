package be.ugent.psb.setas.independent_parsers.Rankings;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class CorrectRankingFiles {

private CommonFunctions cmmFunct;

	
    public ArrayList<Integer> printComplRanks(ArrayList<Integer> nonRepranks, int maximumRank){
	
    	ArrayList<Integer> comp_ranks = new ArrayList<Integer>();
    	
    	for (int i=1; i<maximumRank;i++){
    		
    		if(! cmmFunct.searchIntList_boolean(nonRepranks,i)){
    		comp_ranks.add(i);}
    	}
    	
    	return comp_ranks;
    }
	
	public static void main(String [] args){
		
		CorrectRankingFiles crrRank = new CorrectRankingFiles();
		CommonFunctions cmmFunc = new CommonFunctions();
		ArrayList<Integer> positive_ranks = cmmFunc.read1ColFile_Int("");
	
		//old ones: new one has ~474,000
		int maximumRank_goObo = 250421; //go.obo
		int maximumRank_goBasic = 240698;
		int maximumRank_goSlimPlants = 8726;
		
		ArrayList<Integer> nonRep_ranks = cmmFunc.nonRepetetive_IntList(positive_ranks);
		ArrayList<Integer> comp_ranks = crrRank.printComplRanks(nonRep_ranks, maximumRank_goObo);
		
		for(int i=0; i< nonRep_ranks.size();i++){
			System.out.println(nonRep_ranks.get(i));
		}
	
//		for(int j=0; j<comp_ranks.size();j++){
//			System.out.println(comp_ranks.get(j));
//		}

		
		
		
	}
}
