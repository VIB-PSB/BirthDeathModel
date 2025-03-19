package be.ugent.psb.setas.bdmodel.parsers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class ReadComOutput_optValues {
	
	private ArrayList<String> gfIDs;
	
	public ArrayList<String> getGfIDs() {
		return gfIDs;
	}
	
	public List<List<Double>> readInputFile(String filename) {

		FileReader fin = null;
		try {
			fin = new FileReader(filename);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		/* Skip the first line */
		Scanner sc = new Scanner(fin);
		sc.nextLine();

		List<List<Double>> table = new ArrayList<List<Double>>();
		gfIDs = new ArrayList<String>();
		
		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			String[] myStrings = line.split("\t");

			/* Now we have the numbers as strings in the array "numbers", so: */
					
			List<Double> chunks = new ArrayList<Double>();
			
			// column 0: GF-ID , 1: root size*, 2: lambda*, 3: Lk*
			gfIDs.add(myStrings[0]);
			
			for (int i = 1; i < myStrings.length; i++) {

				double parsed = Double.parseDouble(myStrings[i]);
				chunks.add(parsed);
			}
			table.add(chunks);

		}
		sc.close();
		return table;
	}
	

}
