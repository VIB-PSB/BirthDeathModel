package be.ugent.psb.setas.independent_parsers.Copy;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class CopyRepeatedRecords {
	
	private ArrayList<String> gfIDs_Rep;
	private ArrayList<String> gfIDs_Key;
	
	 public ArrayList<String> getGfIDs_Rep() {
		 return gfIDs_Rep;
		 }
	 
	 public ArrayList<String> getGfIDs_Key() {
		 return gfIDs_Key;
		 }

	public List<List<String>> readMapFile(String mapFile) {

		FileReader fin = null;
		try {
			fin = new FileReader(mapFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		Scanner sc = new Scanner(fin);

		List<List<String>> map = new ArrayList<List<String>>();

		gfIDs_Rep = new ArrayList<String>();
		gfIDs_Key = new ArrayList<String>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();
			
			if(!line.isEmpty()){
			String[] chunks = line.split("\t");

			gfIDs_Rep.add(chunks[0]);
			gfIDs_Key.add(chunks[1]);

			List<String> ls = new ArrayList<String>();
			ls.add(chunks[0]);
			ls.add(chunks[1]);

			map.add(ls);
			}
		}
		sc.close();
		return map;
	}
	
//	public String searchMap(List<List<String>> map, String probe){
//		
//		for(List<String> ls : map){
//			
//			// careful, now the probe must be in the second column: [1]
//			if(ls.get(1).equals(probe)){
//				//and the redult in the first column [0]
//				// if this is not the case, change 0 and 1
//				return ls.get(0);
//			}
//		}
//		
//		return null;
//	}
	
	/* combinedResultFile is the result of the simulations*/
	public List<Double> search_GFID(String combinedResultFile,String key) {
		
		FileReader fin = null;
		try {
			fin = new FileReader(combinedResultFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		/* The first line is a header: */
		Scanner sc = new Scanner(fin);
		sc.nextLine();
		
		List<Double> parameters = new ArrayList<Double>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();
			String[] chunks = line.split("\t");
			
			if(chunks[0].equals(key)){
				
				for (int i = 1; i < chunks.length; i++) {
					double parsed = Double.parseDouble(chunks[i]);	
					parameters.add(parsed);
				}
				sc.close();
				return parameters;
			
			
			}

		}
		sc.close();
		
		return parameters;
		
	}
	
	public static void main(String[] args){
		
		// to creat the map for 9178 core angiosperms with different GF-IDs than 3116 Ortho (Eud or Mon), we wrote a python script 
		
		CopyRepeatedRecords cpRepRec = new CopyRepeatedRecords();
		
		//just run the following line to fill up gfIDs_rep and _key
		List<List<String>> map = cpRepRec.readMapFile("/home/setas/git/StochasticBD/src/Files/GeneFamilyCounts/repeatedGF_ID_Map_1484Ortho_Monocots");
		
		ArrayList<String> gfIDs_Rep = cpRepRec.getGfIDs_Rep();
		ArrayList<String> gfIDs_Key = cpRepRec.getGfIDs_Key();
	
		for(int i =0 ; i< gfIDs_Rep.size(); i++){
			
			// we know that the two lists are corresponding each other at the same index i
			String rep = gfIDs_Rep.get(i);
			String key = gfIDs_Key.get(i);
			
		System.out.print(rep+"\t");
		
		List<Double> GFrecord = cpRepRec.search_GFID("/home/setas/git/StochasticBD/src/Files/GeneFamilyCounts/combinedOutput_1484Ortho_MonocotsMGCF5_CPMpval", key);
		
		for(double d : GFrecord){
			
			System.out.print(d+"\t");
		}
		
		System.out.println();
		
		}
		
	}
}
