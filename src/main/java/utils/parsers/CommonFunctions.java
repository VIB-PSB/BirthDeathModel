package utils.parsers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class CommonFunctions {
	
	public List<List<String>> readMapFile(String mapFile) {

		FileReader fin = null;
		try {
			fin = new FileReader(mapFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);
		List<List<String>> map = new ArrayList<List<String>>();

		while (sc.hasNextLine()) {
			String line = sc.nextLine();

			if (!line.isEmpty()) {
				
				String[] chunks = line.split("\t");
				List<String> ls = new ArrayList<String>();
				
				for(int i=0; i<chunks.length;i++){
				ls.add(chunks[i]);
				}

				map.add(ls);
			}
		}
		sc.close();
		return map;
	}	
}
