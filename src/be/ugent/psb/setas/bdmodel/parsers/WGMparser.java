package be.ugent.psb.setas.bdmodel.parsers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class WGMparser {

	/**
	 * reads a file containing all WGDs on the tree
	 * 
	 * @param filename
	 * @return a list of WGDs to be built via addWGD and InsertWGD in class
	 *         Node.
	 *         
	 *         Warning: If there are multiple WGMs on a branch, the order in which they appear in the text file is important: 
	 *         first the older events and then the younger ones
	 */
	public List<List<String>> readInputFile(String filename) {

		FileReader fr = null;
		try {
			fr = new FileReader(filename);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fr);

		List<List<String>> table = new ArrayList<List<String>>();
		
		while (sc.hasNextLine()) {

			String line = sc.nextLine();
			String[] words = line.split(",");

			List<String> info = new ArrayList<String>();

			for (int i = 0; i < words.length; i++) {
				info.add(words[i]);
			}
			table.add(info);

		}
		sc.close();
		return table;
	}
	
	
	public ArrayList<String> readRetentionRates_String(String filename) {

		ArrayList<String> myStrings = new ArrayList<String>();

		FileReader fin = null;
		try {
			fin = new FileReader(filename);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);

		while (sc.hasNextLine()) {

			String line = sc.nextLine();
			myStrings.add(line);
			
		}
		sc.close();
		return myStrings;
	}

//	 public static void main(String args[]) {
//	 WGM wgm = new WGM();
//	 List<List<String>> list = wgm
//	 .readInputFile("/home/setas/workspace/BirthDeathModel/src/Files/wgd.txt");
//	 System.out.println(list.size());
//	 System.out.println(list.get(0));
//	 System.out.println(list.get(1));
//	
//	 }
}
