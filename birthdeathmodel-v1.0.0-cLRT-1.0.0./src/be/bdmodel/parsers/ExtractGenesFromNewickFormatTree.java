package be.ugent.psb.setas.bdmodel.parsers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class ExtractGenesFromNewickFormatTree {

	public List<List<String>> readMapFile_newickTreesZhen(String mapFile) {

		FileReader fin = null;
		try {
			fin = new FileReader(mapFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);
		sc.nextLine();

		List<List<String>> map = new ArrayList<List<String>>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			if (!line.isEmpty()) {

				String[] chunks = line.split("\t");

				List<String> ls = new ArrayList<String>();

				// we need only 2 columns ..> not anymore!
				ls.add(chunks[0]); //GF-ID in case of break up _1 and _2
				ls.add(chunks[1]); //GF-Id
				ls.add(chunks[2]); //Newick-Tree

				map.add(ls);
			}
		}
		sc.close();
		return map;
	}

	public static void main(String[] args) {

		ExtractGenesFromNewickFormatTree extG = new ExtractGenesFromNewickFormatTree();
		List<List<String>> map = extG.readMapFile_newickTreesZhen("/home/setas/git/StochasticBD/src/Files/genefamily.angiosperm.core.gftree.txt");

		for (List<String> ls : map) {

			String gfID_incaseSplitted = ls.get(0);
			String gfID = ls.get(1);
			String newickTree = ls.get(2);

			String[] nodes = newickTree.split("(\\,)|(\\()|(\\))|(;)"); //node = Atha|AT5G24120; node[1] = AT5G24120

			for (String n : nodes) {

				if (!n.isEmpty()) {
					
					if( n.split("\\|").length > 1){ //if it can not split, the whole thing would be one part and the lenght would be 1
						
					if( gfID_incaseSplitted.split("_").length > 1){	
					System.out.println(gfID_incaseSplitted.split("EE")[1]+"\t"+ n.split("\\|")[1]);
					}
					
					else{
						System.out.println(gfID+"\t"+ n.split("\\|")[1]);
					}
					
					}
				}
			}

		}

	}

}
