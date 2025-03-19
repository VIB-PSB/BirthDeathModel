package be.ugent.psb.setas.bdmodel.parsers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class SeqParser {

	public List<List<String>> readInputFile(String filename) {

		FileReader fin = null;
		try {
			fin = new FileReader(filename);
		} catch (FileNotFoundException e) {
			e.printStackTrace();}
		Scanner sc = new Scanner(fin);
		List<List<String>> table = new ArrayList<List<String>>();

		while (sc.hasNextLine()) {
			String line = sc.nextLine();
			String[] strings = line.split("\t");

			// Now, you have the strings in the array, so:
			List<String> nums = new ArrayList<String>();

			for (int i = 0; i < strings.length; i++) {
				String parsed = strings[i];
				nums.add(parsed);
			}
			table.add(nums);
		}
		sc.close();
		return table;
	}
	
	public List<List<String>> sortTable(List<List<String>> tableOne) {

		List<List<String>> sortedTable = new ArrayList<List<String>>();
		
		boolean[] visited = new boolean[tableOne.size()];
//		System.out.println(visited.length);
		
//		System.out.println(tableOne.size());
		
		for (int i = 0; i < tableOne.size(); i++) {
			
			if (visited[i] == false) {
				
				String probe = tableOne.get(i).get(2);
				
				sortedTable.add(tableOne.get(i));
				
				visited[i] = true;
				
				for (int j = i+1 ; j < tableOne.size()-1; j++) {

                   boolean b= tableOne.get(j).get(2).equals(probe);
					
					if (visited[j] == false 
							&& b==true) {
						
						sortedTable.add(tableOne.get(j));

						visited[j] = true;
					}

				}
			}
		}

		return sortedTable;
	}
	
//	public static void main(String[] args) {
//		// TODO Auto-generated method stub
//		
//		SeqParser sqp = new SeqParser();
//		
//		List<List<String>> table = sqp.readInputFile("/home/setas/workspace/BirthDeathModel/src/Files/sequences.txt");
//
//		
//		List<List<String>> tableTwo = sqp.sortTable(table);
////		
//		for(int i=0; i<tableTwo.size();i++){
//		System.out.println(tableTwo.get(i));}
//	}

}
