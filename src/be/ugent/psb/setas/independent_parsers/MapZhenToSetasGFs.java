package be.ugent.psb.setas.independent_parsers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class MapZhenToSetasGFs {

	private ArrayList<String> gfIDs_zhen;
	private ArrayList<String> gfIDs_setas;

	public ArrayList<String> getGfIDs_zhen() {return gfIDs_zhen;}
	public ArrayList<String> getGfIDs_setas() {return gfIDs_setas;}

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

	public String searchMap(List<List<String>> map, String probe) {
		for (List<String> ls : map) {
			// careful, now the probe must be in the column[0]
			if (ls.get(0).equals(probe)) {
				// and the result in the column [1] if this is not the case, change 0 and 1
				return ls.get(1);
			}
		}
		return null;
	}

	public ArrayList<String> readOneColumnFile(String filename) {

		ArrayList<String> gfIDs = new ArrayList<String>();

		FileReader fin = null;
		try {
			fin = new FileReader(filename);
		} catch (FileNotFoundException e) {
			e.printStackTrace();}
		Scanner sc = new Scanner(fin);
		while (sc.hasNextLine()) {
			String line = sc.nextLine();
			gfIDs.add(line);
		}
		sc.close();
		return gfIDs;
	}

	public static void main(String[] args) {

		MapZhenToSetasGFs mapzs = new MapZhenToSetasGFs();

		ArrayList<String> zhenGFids_inOrderLamda = mapzs
				.readOneColumnFile("/home/setas/git/StochasticBD/src/Files/9178coreGF_inOrderLam");
		List<List<String>> map = mapzs.readMapFile("/home/setas/git/StochasticBD/src/Files/zhen_setas_Mapresults");

//		ArrayList<String> gfIDs_zhen = mapzs.getGfIDs_zhen();
//		ArrayList<String> gfIDs_setas = mapzs.getGfIDs_setas();

		for (String zhenID : zhenGFids_inOrderLamda) {
			String setasID = mapzs.searchMap(map, zhenID);
			if (setasID != null) {
				System.out.println(zhenID + "\t" + setasID);
			}
		}

//	    MapZhenToSetasGFs mapzs = new MapZhenToSetasGFs();
//		
//		ArrayList<String> eudGFids_inOrderLamda = mapzs.readOneColumnFile("/home/setas/git/StochasticBD/src/Files/3116GFIDs_EudPrun_SortLam");
//		
//		List<List<String>> map = mapzs.readMapFile("/home/setas/git/StochasticBD/src/Files/zhen_setas_Mapresults");
//		
//	
//		for(String setasID : eudGFids_inOrderLamda){
//			
//			String zhenID = mapzs.searchMap(map, setasID);
//			
//			if(zhenID != null){
//				//we still print setasID becasue we want to have both files in one terminology: ORTHSTMxxx to be able to compare their rankings
//			System.out.println(setasID);}
//			
//		}
	}

}