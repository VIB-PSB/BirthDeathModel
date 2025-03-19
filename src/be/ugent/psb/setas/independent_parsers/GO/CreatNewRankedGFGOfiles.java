package be.ugent.psb.setas.independent_parsers.GO;


import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class CreatNewRankedGFGOfiles {
	
	
	public static void main(String[] args) {

		CommonFunctions cmmFunct = new CommonFunctions();
		
//		String path_rank_GFid ="/home/setas/Desktop/setas/Project1/Results/CompareRankings/PCCranks-Matlab/newRanks.txt";
		String path_rank_GFid ="/home/setas/Desktop/setas/Project1/Results/CompareRankings/RankCorrelations/PPCranks-Matlab/old/newRanks.txt";
		
		ArrayList<String> GFids = cmmFunct.readColX_String(path_rank_GFid, 0);
		int colNumberWithRanks = 8; //Angio-Tau-Mon2-Musa3
		ArrayList<Double> ranks = cmmFunct.readColX_double(path_rank_GFid, colNumberWithRanks); 
		
		String path_GFid_GOid_Desc="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/MannWhitneyUtest/ranks_gfID_goHierarch_Desc_inOrderLam";
//		String path_GFid_GOid_Desc ="/home/setas/git/IndependentParsers/Independent Parsers/src/files/37AngioSperms/rankedFiles/ranks_GFidOrderLam_GOid_Desc_GoSlimPlants_37spe";
		List<List<String>> mapGFidGOidDesc = cmmFunct.readMapFile(path_GFid_GOid_Desc);
		
		for(int i=0; i<GFids.size();i++){
			
			
			String gfIDcurrent = GFids.get(i);
			double rankcurrent = ranks.get(i);
			
//			System.out.println(gfIDcurrent+"\t"+rankcurrent);
			
			for(int j=0; j< mapGFidGOidDesc.size();j++){
				
				List<String> line = mapGFidGOidDesc.get(j);
				
//				System.out.println("line.get(1)  "+line.get(1));
				
				if(line.get(1).equals(gfIDcurrent)){ //0==old-rank and line.get(1)= gfID
					
					System.out.print(rankcurrent+"\t"+gfIDcurrent+"\t");
					
					for(int k=2; k<line.size();k++){ // forget about the previous rank =0 and gf id is already printed =1
						
						System.out.print(line.get(k)+"\t");
					}
					System.out.println();
				}	
			}
		}
		
		
		
		
		
		
	}

}
