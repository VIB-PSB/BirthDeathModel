package be.ugent.psb.setas.independent_parsers.Expression;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class ReadLoessFitExpDataOutputFiles {
	
	CommonFunctions cmf = new CommonFunctions();
	
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
				
				String[] chunks = line.split(",");

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
	
	
	public String[] readCSoutout_returnInOneList(String myFile, boolean header){
		
		List<List<String>> lst = readMapFile(myFile);
		String [] lst_array = new String [lst.size()*3];
		
//		System.out.println(lst_array.length);
		
		
		List<String> firstCol = cmf.readColX_String_Delimiter(myFile, 0, ",", header);
		List<String> secondCol = cmf.readColX_String_Delimiter(myFile, 1, ",", header);
		List<String> thirdCol = cmf.readColX_String_Delimiter(myFile, 2, ",", header);
		
		int j=0;
		for (int i = 0; i < lst.size(); i++) {
				lst_array[j++] = firstCol.get(i);
				lst_array[j++] = secondCol.get(i);
				lst_array[j++] = thirdCol.get(i);
		}
		
		return lst_array;
		
		
	}
	
	public static void main(String [] args){
		
		CommonFunctions cmf2 = new CommonFunctions();
		ReadLoessFitExpDataOutputFiles readLoess = new ReadLoessFitExpDataOutputFiles();
		
		String path_realValue= "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/NewFits/realValues.txt";
		String path_fittedValues="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/NewFits/fitted.txt";
		String path_residues="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/NewFits/residuals.txt";
		
		String path_GF_KsLessThan5 ="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/NewFits/ExpKs_GFKsLess5";
		
		List<List<String>> myMainInfo = cmf2.readMapFile(path_GF_KsLessThan5);
		
		String [] reals = readLoess.readCSoutout_returnInOneList(path_realValue, false);
		String [] fitteds = readLoess.readCSoutout_returnInOneList(path_fittedValues, false);
		String [] residues = readLoess.readCSoutout_returnInOneList(path_residues, false);
		
		
		for(int i=0; i<myMainInfo.size();i++){
			
			List<String> lsCurrent = myMainInfo.get(i);
			
			for(int j=0; j<lsCurrent.size();j++){
				
				System.out.print(lsCurrent.get(j)+"\t");
			}
			
			System.out.println(reals[i]+"\t"+fitteds[i]+"\t"+residues[i]);
		}
		
		
		
		
		
		
		
	}

}
