package be.ugent.psb.setas.independent_parsers.GO;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;
import be.ugent.psb.setas.independent_parsers.Rankings.CompareRankings;

public class PrintAthaGO_ForCombinedOutput_refined {

	public List<List<String>> readMapFile(String mapFile) {

		FileReader fin = null;
		try {
			fin = new FileReader(mapFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);
		List<List<String>> map = new ArrayList<List<String>>();

		ArrayList<String> gfIDs = new ArrayList<String>();
		ArrayList<String> goIDs = new ArrayList<String>();
		ArrayList<String> goDesc = new ArrayList<String>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			if (!line.isEmpty()) {
				String[] chunks = line.split("\t");

				gfIDs.add(chunks[0]);
				goIDs.add(chunks[1]);
				goDesc.add(chunks[2]);

				List<String> ls = new ArrayList<String>();
				ls.add(chunks[0]);
				ls.add(chunks[1]);
				ls.add(chunks[2]);

				map.add(ls);
			}
		}
		sc.close();
		return map;
	}

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

//		PrintAthaGenesGFsInOrder prAthGen = new PrintAthaGenesGFsInOrder();
		CommonFunctions cmmFunct = new CommonFunctions();
		CompareRankings cmpRank = new CompareRankings();

		// ArrayList<String> GF_ids_inOrder = new ArrayList<String>();
		// GF_ids_inOrder =
		// cmmFunct.read1ColFile_String("/home/setas/git/IndependentParsers/Independent Parsers/src/files/9178coreGF_inOrderLam");
		// List<List<String>> map
		// =prAthGen.readMapFile("/home/setas/git/StochasticBD/src/Files/Atha.core.txt");
		// List<List<String>> map = prAthGen
		// .readMapFile("/home/setas/git/IndependentParsers/Independent Parsers/src/files/GoHierarchy_GoObo_Desc");
		// for (String gfID : GF_ids_inOrder) {
		// List<List<String>> found = prAthGen.searchMap(mapGFAtha, gfID);
		// }

		String comOutput_allMeasures = "/home/setas/Desktop/setas/Project1/Results/CompareRankings/combinedOutput_reGen.txt";
		List<List<String>> mapAllMeasures = cmmFunct.readMapFile(comOutput_allMeasures);

		String path_GFAthaMap = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/Tair/addedInforFromTair/GFid_Ath_geneDescrTair_InOrder";
		List<List<String>> mapGFAtha = cmpRank.readMapFile(path_GFAthaMap);
		// List<String> GFids_Atha = cmmFunct.readColX_String(path_GFAthaMap,
		// 0);

		String path_GOids = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/Tair/addedInforFromTair/GFid_Ath_ATHGOSLIM.txt";
		List<List<String>> mapGFAthaGO = cmpRank.readMapFile(path_GOids);

		/** for final combined output file */

		for (int i = 0; i < mapAllMeasures.size(); i++) {

			List<String> gFrecord_allMeasures = mapAllMeasures.get(i);
			String gfID_current = gFrecord_allMeasures.get(0);
			
			System.out.print(gfID_current + "\t");

			List<List<String>> found_AthaDesc_1GF = cmpRank.searchMap(mapGFAtha, gfID_current);
			List<String> saveAthaGenes = new ArrayList<String>();

			List<List<String>> found_AthaGO = cmpRank.searchMap(mapGFAthaGO, gfID_current);
			List<String> goIDs = new ArrayList<String>();
			
			for (int j = 0; j < found_AthaGO.size(); j++) {

				if (found_AthaGO.get(j).size() > 6) {

					String goID = found_AthaGO.get(j).get(6).split(":")[1];

					if (!(goIDs.contains(goID))) {

						goIDs.add(goID);
					}
				}
			}

			/** print the A.thaliana gene IDs */
			if (found_AthaDesc_1GF.size() > 0) {

				for (int k = 0; k < found_AthaDesc_1GF.size(); k++) {

					List<String> line_AthaDesc = found_AthaDesc_1GF.get(k);

					if (!(saveAthaGenes.contains(line_AthaDesc.get(1)))) { /** if it is not repetetive Atha genes,
						// because sometime the same Atha is in two lines,due to source of description */

						saveAthaGenes.add(line_AthaDesc.get(1));
						System.out.print(line_AthaDesc.get(1));

						if (k != found_AthaDesc_1GF.size() - 1) {
							System.out.print(",");
						}

					}

				}

				System.out.print("\t");		
				
				/** print A.thaliana GO IDs*/
				
				if (goIDs.size() == 0) {
					System.out.print("N/A");
				}

				else{
				for (int m = 0; m < goIDs.size(); m++) {

						System.out.print(goIDs.get(m));

						if (m != goIDs.size() - 1) {
							System.out.print(",");
						}
					
				}
				}

				/** print the description, only from the first line corrsponding to first Atha gene. Because other genes have more or less same function*/
				System.out.print("\t"+ found_AthaDesc_1GF.get(0).get(3).split(";")[0]+ "\t");

			}

			else {

				System.out.print("NA" + "\t" + "NA"+"\t"+"NA" +"\t");

			}
			
			for (int j = 1; j < gFrecord_allMeasures.size()-1; j++) {
				System.out.print(gFrecord_allMeasures.get(j) + "\t");
			}

			System.out.print(gFrecord_allMeasures.get(gFrecord_allMeasures.size()-1));
			System.out.println();

		}

	}
}
