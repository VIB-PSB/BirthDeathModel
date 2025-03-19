package be.ugent.psb.setas.bdmodel.parsers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class prunMGFgeneCounts {

	public int searchListString(String prob, ArrayList<String> myList) {

		int index = -1;
		
		for(int i=0; i<myList.size();i++){
			
			if(myList.get(i).equals(prob)){
				
				index = i;
			}
		}
		
		return index;
		
	}
	
	public ArrayList<String> readOneColumnFile(String filename) {

		ArrayList<String> gfIDs = new ArrayList<String>();

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
			gfIDs.add(line.trim());
			
		}
		sc.close();
		return gfIDs;
	}

	public static void main(String[] args) {

		ReadGFcountsFile rgf = new ReadGFcountsFile();
		List<List<Integer>> counts = rgf.read_all(args[0]);
		ArrayList<String> gfIDs = rgf.getGfIDs();
		
		prunMGFgeneCounts pMGF = new prunMGFgeneCounts();
		
		ArrayList<String> MGFids = pMGF.readOneColumnFile(args[1]);
		
		for(String mgf : MGFids){
			
			int index = pMGF.searchListString(mgf, gfIDs);	
			List<Integer> mgfCounts = counts.get(index);
			
			System.out.println(mgf);
			
			for(int mgfC : mgfCounts){
				
				System.out.println(mgfC);
			}
		}
		
		

	}
}
