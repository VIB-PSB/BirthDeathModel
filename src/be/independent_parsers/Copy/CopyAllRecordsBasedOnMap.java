package be.ugent.psb.setas.independent_parsers.Copy;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class CopyAllRecordsBasedOnMap {

	// function 1: read map file:
	public List<List<String>> readMapFile(String mapFile) {

		FileReader fin = null;
		try {
			fin = new FileReader(mapFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);

		List<List<String>> map = new ArrayList<List<String>>();

		ArrayList<String> gfIDs_zhen = new ArrayList<String>();
		ArrayList<String> gfIDs_setas = new ArrayList<String>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			if (!line.isEmpty()) {
				String[] chunks = line.split("\t");

				gfIDs_zhen.add(chunks[0]);
				gfIDs_setas.add(chunks[1]);

				List<String> ls = new ArrayList<String>();
				ls.add(chunks[0]);
				ls.add(chunks[1]);

				map.add(ls);
			}
		}
		sc.close();
		return map;
	}

	// // function 2: serach map for a key GF_ID: return a corresponding GF-ID
	// in
	// // the second file
	// public String searchMap(List<List<String>> map, String probe) {
	//
	// for (List<String> ls : map) {
	//
	// // careful, now the probe must be in the second column [1]
	// if (ls.get(1).equals(probe)) {
	//
	// /*
	// * and the result in the first column [0] if this is not the
	// * case, change 0 and 1
	// */
	// return ls.get(0);
	//
	// }
	// }
	//
	// return null;
	// }

	// function 3: print out all records for a GF_ID from a combined reults
	// file:

	/* combinedResultFile is the result of the simulations */
	public void searchGFid_printRec(String combinedResultFile, String gfID) {

		FileReader fin = null;
		try {
			fin = new FileReader(combinedResultFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);

		while (sc.hasNextLine()) {

			List<Double> parameters = new ArrayList<Double>();
			String line = sc.nextLine();
			String[] chunks = line.split("\t");

			if (chunks[0].equals(gfID)) {

				System.out.print(chunks[0] + "\t");

//				for (int i = 1; i < 4; i++) {
					 for (int i = 1; i < chunks.length; i++) {
					double parsed = Double.parseDouble(chunks[i]);

					System.out.print(parsed + "\t");
					parameters.add(parsed);
				}
			}

		}
		sc.close();
	}

	public void searchGFid_printGeneCounts_all(String geneCountsFile,
			String gfID) {

		FileReader fin = null;
		try {
			fin = new FileReader(geneCountsFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);

		while (sc.hasNextLine()) {

			List<Integer> geneCounts = new ArrayList<Integer>();
			String line = sc.nextLine();
			String[] chunks = line.split("\t");

			if (chunks[0].equals(gfID)) {

				for (int i = 3; i < chunks.length; i++) {
					int parsed = Integer.parseInt(chunks[i]);

					System.out.print(parsed + "\t");
					geneCounts.add(parsed);
				}
			}

		}
		sc.close();
	}
	
	public void read_prntParam_prntGenCounts(String combinedResultFile, String geneCountsFile) {

		FileReader fin = null;
		try {
			fin = new FileReader(combinedResultFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);

		while (sc.hasNextLine()) {

			String line = sc.nextLine();
			String[] chunks = line.split("\t");

			System.out.print(chunks[0] + "\t");

			for (int i = 1; i < chunks.length; i++) {
				double parsed = Double.parseDouble(chunks[i]);
				System.out.print(parsed + "\t");
			}
			
			searchGFid_printGeneCounts_all(geneCountsFile, chunks[0]);
			System.out.println();

		}
		sc.close();
	}

	public void searchGFid_printGeneCounts_selective(String geneCountsFile,
			String gfID, int[] selectedIndex) {

		FileReader fin = null;
		try {
			fin = new FileReader(geneCountsFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);

		while (sc.hasNextLine()) {

			List<Integer> geneCounts = new ArrayList<Integer>();
			String line = sc.nextLine();
			String[] chunks = line.split("\t");

			if (chunks[0].equals(gfID)) {

				System.out.print(chunks[0] + "\t");

				for (int i = 3; i < chunks.length; i++) {
					int parsed = Integer.parseInt(chunks[i]);
					geneCounts.add(parsed);
				}

				for (int index : selectedIndex) {

					System.out.print(geneCounts.get(index)+"\t");
				}
			}

		}
		sc.close();
	}

	public static void main(String[] args) {

//		CopyAllRecordsBasedOnMap cpAllRec = new CopyAllRecordsBasedOnMap();
//		List<List<String>> map = cpAllRec
//				.readMapFile("/home/setas/git/StochasticBD/src/Files/zhenInOrder_Tosetas_notNull");
//
//		int[] selectedIndex = new int[10];
//
//		selectedIndex[0] = 0;
//		selectedIndex[1] = 4;
//		selectedIndex[2] = 5;
//		selectedIndex[3] = 11;
//		selectedIndex[4] = 13;
//		selectedIndex[5] = 17;
//		selectedIndex[6] = 18;
//		selectedIndex[7] = 23;
//		selectedIndex[8] = 25;
//		selectedIndex[9] = 26;

//		for (List<String> ls : map) {
//
//			// ls.get(0) = zhe's ids, ls.get(1) = stmae-ids
//			if (ls.get(0) != null && ls.get(1) != null) {
//
//				// cpAllRec.searchGFid_printRec("/home/setas/git/StochasticBD/src/Files/comSortLam_37speMGCF5_CPMpval_9178coreGF",
//				// ls.get(0));
//				// cpAllRec.searchGFid_printRec("/home/setas/git/StochasticBD/src/Files/combined_3116OrthGF_CPM_Pvals_SortLam.txt",
//				// ls.get(1));
//
//				cpAllRec.searchGFid_printGeneCounts_selective(
//						"/home/setas/git/StochasticBD/src/Files/GeneFamilyCounts/37spe-MGCF5-9178core-OrderNewickTree",
//						ls.get(0), selectedIndex);
//				System.out.println();
//
//				cpAllRec.searchGFid_printGeneCounts_all(
//						"/home/setas/git/StochasticBD/src/Files/GeneFamilyCounts/OrthoGF_inAll14spe_EudictosMGCF5.txt",
//						ls.get(1));
//				System.out.println();
//			}
//
//			System.out.println();
//
//		}
		
		CopyAllRecordsBasedOnMap cpAllRec = new CopyAllRecordsBasedOnMap();
		
		String combinedResultFile = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/First709GF_CombOutput";
		String geneCountsFile ="/home/setas/git/StochasticBD/src/Files/GeneFamilyCounts/37spe-MGCF5-9178core-OrderNewickTree";
		cpAllRec.read_prntParam_prntGenCounts(combinedResultFile, geneCountsFile);
		

	}

}
