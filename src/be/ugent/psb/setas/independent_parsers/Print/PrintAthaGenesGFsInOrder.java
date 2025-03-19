package be.ugent.psb.setas.independent_parsers.Print;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;
import be.ugent.psb.setas.independent_parsers.Rankings.CompareRankings;

public class PrintAthaGenesGFsInOrder {

	public List<List<String>> searchMap(List<List<String>> map, String probe) {

		List<List<String>> found = new ArrayList<List<String>>();

		for (List<String> ls : map) {
			// careful, now the probe must be in the column: [0]
			if (ls.get(0).equals(probe)) {

				found.add(ls);

			}
		}

		return found;

	}

	public static void main(String[] args) {

		CommonFunctions cmmFunct = new CommonFunctions();
		CompareRankings cmpRank = new CompareRankings();

		String comOutput_allMeasures = "/home/setas/Desktop/setas/Project1/Results/CompareRankings/comOut_AllMeasures.txt";
		List<List<String>> mapAllMeasures = cmmFunct
				.readMapFile(comOutput_allMeasures);

		String path_GFAthaMap = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/Tair/addedInforFromTair/GFid_Ath_geneDescrTair_InOrder";
		List<List<String>> mapGFAtha = cmpRank.readMapFile(path_GFAthaMap);

		String path_GOids = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/Tair/addedInforFromTair/GFid_Ath_ATHGOSLIM.txt";
		List<List<String>> mapGFAthaGO = cmpRank.readMapFile(path_GOids);

		for (int i = 0; i < mapAllMeasures.size(); i++) {

			List<String> gfRecord_allMeasures = mapAllMeasures.get(i);
			String gfID_current = gfRecord_allMeasures.get(0);
			
			System.out.print(gfID_current+"\t");

			List<List<String>> found_AthaDesc = cmpRank.searchMap(mapGFAtha,gfID_current);
			List<String> saveNonRepAthaGenes = new ArrayList<String>();

			List<List<String>> found_AthaGO = cmpRank.searchMap(mapGFAthaGO,gfID_current);
			List<String> saveNonRepGOIDs = new ArrayList<String>();
			
			String description="NA";

			for (int j = 0; j < found_AthaGO.size(); j++) {

				if (found_AthaGO.get(j).size() > 6) {

					String goID = found_AthaGO.get(j).get(6).split(":")[1];

					if (!(saveNonRepGOIDs.contains(goID))) {

						saveNonRepGOIDs.add(goID);
					}
				}
			}

			if (found_AthaDesc.size() > 0) {
				
			for (int j = 0; j < found_AthaDesc.size(); j++) {

					String AthaGeneID = found_AthaDesc.get(j).get(1);

					if (!(saveNonRepAthaGenes.contains(AthaGeneID ))) {

						saveNonRepAthaGenes.add(AthaGeneID );
					}
				}
			
			description = found_AthaDesc.get(0).get(3);
			}	
			
			if(saveNonRepAthaGenes.size() >0){
			for(String athID : saveNonRepAthaGenes){
				System.out.print(athID+",");
			}
			}
			
			else{
				System.out.print("NA");
			}
			
			System.out.print("\t"+description+"\t");
			
			if(saveNonRepGOIDs.size()> 0){

			for(String go : saveNonRepGOIDs){
				
				System.out.print(go+",");
			}
			}
			
			else{
				System.out.print("NA");
			}
			
			System.out.println();
			
		}
	}
}
