package be.ugent.psb.setas.bdmodel.parsers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class ReadGenePairDistances {

	/* To read allAnchorPoint file and save all of the numbers in a table */
	public List<List<Double>> readInputFile(String filename) {
	FileReader fin = null;
	try {
		fin = new FileReader(filename);
	} catch (FileNotFoundException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}

	/* The first line is a header: */
	Scanner sc = new Scanner(fin);

    List<List<Double>> table = new ArrayList<List<Double>>();

	while (sc.hasNextLine()) {

		String line = sc.nextLine();
		String[] numbers = line.split("\t");
		
		List<Double> nums = new ArrayList<Double>();
		
		/* The first two columns are gene pairs: */
		for (int i = 2; i < numbers.length; i++) {
			
		/* Now we have the numbers as strings in the array "numbers", so: */
			
			
		if(! numbers[i].equalsIgnoreCase("")){
		double parsed = Double.parseDouble(numbers[i]);
			
		nums.add(parsed);}
		
		}
		/* To avoid pairs with not enough info for t and Ks and else! */
//		if(nums.size() > 6){
		table.add(nums);
//		}
	}
	sc.close();
	return table;
}
	/* To filter the anchorpoints pairs with big Ks */
	public List<List<Double>> filterKsgreaterThan(List<List<Double>> firstTable , double KsThreshold) {
		
		List<List<Double>> filteredLines = new ArrayList<List<Double>>();
		
		for(int i=0; i< firstTable.size(); i++){
			
			if(firstTable.get(i).get(5)< KsThreshold){
				
				filteredLines.add(firstTable.get(i));
			}
			
		}
		
		return filteredLines;
		
	}
	
	/* To filter the anchorpoints pairs out of a specific range for Ks */
	public List<List<Double>> filterKsOutoFRange(List<List<Double>> firstTable , double KsMin, double KsMax) {
		
		List<List<Double>> filteredLines = new ArrayList<List<Double>>();
		
		for(int i=0; i< firstTable.size(); i++){
			
			if(firstTable.get(i).get(5)< KsMax && firstTable.get(i).get(5) > KsMin){
				
				filteredLines.add(firstTable.get(i));
			}
			
		}
		
		return filteredLines;
		
	}

	public static void main(String[] args) {

		List<List<Double>> table = new ArrayList<List<Double>>();
		List<List<Double>> fTable = new ArrayList<List<Double>>();
		
		ReadGenePairDistances rgpd = new ReadGenePairDistances();

		table = rgpd.readInputFile("/home/setas/workspace/BirthDeathModel/src/Files/GossypiumAllDatedAnchorPointsTab");
		
		fTable = rgpd.filterKsgreaterThan(table, 5);
//		fTable = rgpd.filterKsOutoFRange(table,0.1,0.3);
		
		NumberFormat format	= NumberFormat.getInstance();
		format.setMaximumFractionDigits(5);
		
//		double sum =0;
		
		for(int i=0; i< fTable.size(); i++){
			
			String stringFormNumber = format.format(fTable.get(i).get(0)/6);
			double d = Double.parseDouble(stringFormNumber);
			
//			sum += fTable.get(i).get(0)/6;
			System.out.print(d+";");
		}
		
//		double avg = (sum*1.000)/(fTable.size()*1.000);
//		System.out.println(avg);
//		
//		double sumSquares=0;
//		for(int i=0; i< fTable.size(); i++){
//			
//			sumSquares+= Math.pow((fTable.get(i).get(0)/6)- avg ,2);
//		}
//	
//		double variance = sumSquares/fTable.size();
//		double std = Math.sqrt(variance);
//		
//		System.out.println(std);
	}
	
}
