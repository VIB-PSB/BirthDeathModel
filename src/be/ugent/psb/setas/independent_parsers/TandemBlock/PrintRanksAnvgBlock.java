package be.ugent.psb.setas.independent_parsers.TandemBlock;

import java.util.ArrayList;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class PrintRanksAnvgBlock {
	

	
	
	public static void main(String[] args) {

		CommonFunctions cmmFunc = new CommonFunctions();
		
//		String myPath="/home/setas/Desktop/AllGFs-All37spe-AvgBlockPercentOrdered.txt";
		
		String myPath="/home/setas/Desktop/setas/Project1/Results/CompareRankings/AvgDiffIdealProfiles/avgDistToIdealProfsAll37SpeOrdered.txt";
		
		ArrayList<String> gfIDs = cmmFunc.readColX_String(myPath, 0);
		ArrayList<Double> avgBlockPercent = cmmFunc.readColX_double(myPath, 1);
//		ArrayList<Double> avgDistance = cmmFunc.readColX_double(myPath, 1);
		
		
		
		
		int rank=1;
		int goDownInRank =0;
		
		System.out.println(rank+"\t"+gfIDs.get(0)+"\t"+avgBlockPercent.get(0));
		
		for(int i=1; i<gfIDs.size();i++){
			
			if(Math.abs(avgBlockPercent.get(i)-avgBlockPercent.get(i-1)) < 0.00001){
				
//				System.out.println(avgBlockPercent.get(i)-avgBlockPercent.get(i-1) );
				
			goDownInRank +=1;
			
			System.out.println(rank+"\t"+gfIDs.get(i)+"\t"+avgBlockPercent.get(i));
			}
			
			else{
				
				
//				System.out.println("entered else");
				goDownInRank +=1;
				rank = 1+ goDownInRank;
				
				System.out.println(rank+"\t"+gfIDs.get(i)+"\t"+avgBlockPercent.get(i));
			}
			
			

		}
		
		
	}

}
