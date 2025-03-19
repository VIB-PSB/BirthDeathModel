package be.ugent.psb.setas.bdmodel.parsers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class ReadGFcountsFile {

	private ArrayList<String> gfIDs;
	private ArrayList<String> gfIDs_unique;
	private ArrayList<String> gfIDs_unique2;
	private ArrayList<String> gfIDs_repeated;

	/**
	 * reads gene family counts at the leaves from an input file
	 * 
	 * @param filename
	 * @return list of integers corresponding to gene family counts at the
	 *         leaves, only for the non-repetetive gene counts
	 * 
	 *         Warning: The gene count file must have species in the order of
	 *         pre-order traverse of the tree
	 */

	public ArrayList<String> getGfIDs_unique() {
		return gfIDs_unique;
	}

	public ArrayList<String> getGfIDs() {
		return gfIDs;
	}

	public List<List<Integer>> read_all(String filename) {

		FileReader fin = null;
		try {
			fin = new FileReader(filename);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		/* The first line is a header: */
		Scanner sc = new Scanner(fin);
		sc.nextLine();

		List<List<Integer>> gfCounts = new ArrayList<List<Integer>>();

		gfIDs = new ArrayList<String>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();
			String[] chunks = line.split("\t");

			List<Integer> geneCounts = new ArrayList<Integer>();

			// column 0: GF-ID , 1: num of genes, 2: num of Spe
			for (int i = 3; i < chunks.length; i++) {
				int parsed = Integer.parseInt(chunks[i]);
				geneCounts.add(parsed);
			}
			gfCounts.add(geneCounts);
			gfIDs.add(chunks[0]);
		}
		sc.close();
		return gfCounts;
	}

	public List<List<Integer>> read_unique(String filename) {

		FileReader fin = null;
		try {
			fin = new FileReader(filename);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		/* The first line is a header: */
		Scanner sc = new Scanner(fin);
		sc.nextLine();

		List<List<Integer>> uniqueCounts = new ArrayList<List<Integer>>();

		gfIDs_unique = new ArrayList<String>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();
			String[] chunks = line.split("\t");

			List<Integer> geneCounts = new ArrayList<Integer>();

			// column 0: GF-ID , 1: num of genes, 2: num of Sp

			for (int i = 3; i < chunks.length; i++) {
				int parsed = Integer.parseInt(chunks[i]);
				geneCounts.add(parsed);
			}

			if (!uniqueCounts.contains(geneCounts)) {

				uniqueCounts.add(geneCounts);
				gfIDs_unique.add(chunks[0]);
			}

		}
		sc.close();
		return uniqueCounts;
	}

	/* creates a map file to be used in coperepeatedRecords */
	public List<List<Integer>> read_Repetetive(String filename) {

		FileReader fin = null;
		try {
			fin = new FileReader(filename);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		/* The first line is a header: */
		Scanner sc = new Scanner(fin);
		sc.nextLine();

		List<List<Integer>> uniqueCounts2 = new ArrayList<List<Integer>>();
		List<List<Integer>> repeatedCounts = new ArrayList<List<Integer>>();

		gfIDs_unique2 = new ArrayList<String>();
		gfIDs_repeated = new ArrayList<String>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();
			String[] chunks = line.split("\t");

			List<Integer> geneCounts = new ArrayList<Integer>();
			for (int i = 3; i < chunks.length; i++) {
				int parsed = Integer.parseInt(chunks[i]);
				geneCounts.add(parsed);
			}

			if (!uniqueCounts2.contains(geneCounts)) {

				uniqueCounts2.add(geneCounts);
				gfIDs_unique2.add(chunks[0]);
			}

			else { /* These are GFs with repeated gene counts */

				int keyIndex = uniqueCounts2.indexOf(geneCounts);
				String keyGFID = gfIDs_unique2.get(keyIndex);

				repeatedCounts.add(geneCounts);
				gfIDs_repeated.add(chunks[0]);

				System.out.println(chunks[0] + "\t" + keyGFID);
				System.out.println();
			}

		}
		sc.close();
		return repeatedCounts;
	}

	/* To use in RemoveColumnsWithZero */
	public List<List<Integer>> reversetable(List<List<Integer>> tab) {

		int length = tab.get(0).size();
		// System.out.println("length: "+length);

		List<List<Integer>> tableReverse = new ArrayList<List<Integer>>();
		for (int j = 0; j < length; j++) {
			tableReverse.add(new ArrayList<Integer>());
		}

		for (int i = 0; i < tab.size(); i++) {
			for (int j = 0; j < tab.get(i).size(); j++) {
				List<Integer> row = tableReverse.get(j);
				row.add(tab.get(i).get(j));
			}
		}
		return tableReverse;
	}

	// public List<List<Integer>> uniqueRows(List<List<Integer>> table){
	//
	// List<List<Integer>> result = new ArrayList<List<Integer>>();
	//
	// for (List<Integer> row : table){
	// if (!result.contains(row)){
	// result.add(row);
	// }
	// }
	//
	// return result;
	// }

	public List<List<Integer>> removeColumnsWithZero(List<List<Integer>> tb) {

		List<List<Integer>> tbReverse = reversetable(tb);

		int nomOfRows = tbReverse.size();
		int lengtOfRows = tbReverse.get(0).size();

		List<List<Integer>> tb2 = new ArrayList<List<Integer>>();

		int[] containZero = new int[nomOfRows];

		for (int i = 0; i < nomOfRows; i++) {

			for (int j = 0; j < lengtOfRows; j++) {

				if (tbReverse.get(i).get(j) == 0) {

					// System.out.println("this column has Zeros: " + i
					// + " hence corresponding species should be ignored");
					containZero[i] = -1;
					break;
				}

			}
		}

		for (int k = 0; k < nomOfRows; k++) {

			if (containZero[k] != -1) {
				tb2.add(tbReverse.get(k));
			}
		}

		return reversetable(tb2);
	}

	public int[][] matrixOfGFdata(List<List<Integer>> tb) {

		int[][] matrixOfGFdata = new int[tb.size()][tb.get(0).size()];

		for (int gf = 0; gf < tb.size(); gf++) {

			List<Integer> li = tb.get(gf);

			for (int m = 0; m < li.size(); m++) {

				matrixOfGFdata[gf][m] = li.get(m);
			}

		}

		return matrixOfGFdata;
	}

//	public static void main(String[] args) {
//
//		ReadGFcountsFile rgf = new ReadGFcountsFile();
//
//		List<List<Integer>> uniqueCounts = rgf
//				.read_unique("/home/setas/git/StochasticBD/src/Files/GeneFamilyCounts/3116Ortho_all14_inOrderMGCF5pruned14.txt");
//
//		List<List<Integer>> repeatedCounts = rgf
//				.read_Repetetive("/home/setas/git/StochasticBD/src/Files/GeneFamilyCounts/3116Ortho_all14_inOrderMGCF5pruned14.txt");
//
//	}

}
