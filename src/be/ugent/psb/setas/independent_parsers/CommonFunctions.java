package be.ugent.psb.setas.independent_parsers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class CommonFunctions {
	
	public ArrayList<String> read1ColFile_String(String filename) {

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
	
	public ArrayList<Integer> read1ColFile_Int(String filename) {

		ArrayList<Integer> myInts = new ArrayList<Integer>();

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
			myInts.add(Integer.parseInt(line));
			
		}
		sc.close();
		return myInts;
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
	
	
	
//	public List<List<String>> readAnnoFile(String mapFile) {
//
//		FileReader fin = null;
//		try {
//			fin = new FileReader(mapFile);
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
//
//		Scanner sc = new Scanner(fin);
//
//		List<List<String>> map = new ArrayList<List<String>>();
//
//		while (sc.hasNextLine()) {
//
//			String line = sc.nextLine();
//
//			if (!line.isEmpty()) {
//				
//				String[] chunks = line.split("=");
//
//				List<String> ls = new ArrayList<String>();
//				
//				for(int i=0; i<chunks.length;i++){
//				ls.add(chunks[i]);
//				}
//
//				map.add(ls);
//			}
//		}
//		sc.close();
//		return map;
//	}
	
	/* The specified column (counitng from 0) should contain doubles or ints*/
	public List<List<String>> subMap(String mapFile, int colNum, double upperLimit) {

		FileReader fin = null;
		try {
			fin = new FileReader(mapFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);

		List<List<String>> sub_map = new ArrayList<List<String>>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			if (!line.isEmpty()) {
				
				String[] chunks = line.split("\t");
				List<String> ls = new ArrayList<String>();
				
				for(int i=0; i<chunks.length;i++){
				ls.add(chunks[i]);
				}
				
				double colContent = Double.parseDouble(chunks[colNum]);
                
				if(colContent < upperLimit){
				sub_map.add(ls);}
			}
		}
		sc.close();
		return sub_map;
	}
	
	public ArrayList<Integer> nonRepetetive_IntList(ArrayList<Integer> myArrayListOfIntegers){
		
		ArrayList<Integer> nonRep = new ArrayList<Integer>();
		
		int currentInt= myArrayListOfIntegers.get(0);
		
		nonRep.add(currentInt);
		
		for(int i=1; i<myArrayListOfIntegers.size(); i++){
			
			int newInt= myArrayListOfIntegers.get(i);
			
			if(newInt != currentInt){
				
				nonRep.add(newInt);
				currentInt = newInt;
			}
		}
		
		return nonRep;
	}
		
	public int[] convertListToArray(ArrayList<Integer> ls) {

		int[] arr = new int[ls.size()];

		for (int i = 0; i < ls.size(); i++) {

			arr[i] = ls.get(i).intValue();
		}
		return arr;
	}
	public double[] convertListToArray_double(ArrayList<Double> ls) {

		double[] arr = new double[ls.size()];

		for (int i = 0; i < ls.size(); i++) {

			arr[i] = ls.get(i).doubleValue();
		}
		return arr;
	}
	
	public ArrayList<String> readColX_String(String fileName, int colNumber) {

		FileReader fin = null;
		try {
			fin = new FileReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);
		ArrayList<String> oneColumn = new ArrayList<String>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			if (!line.isEmpty()) {

				String[] chunks = line.trim().split("\t");
				oneColumn.add(chunks[colNumber]);
			}
		}
		sc.close();
		return oneColumn;
	}
	

	
	public ArrayList<String> readColX_String_Delimiter(String fileName, int colNumber, String delimiter, boolean header) {

		FileReader fin = null;
		try {
			fin = new FileReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);
		if(header){
			sc.nextLine();
		}
		ArrayList<String> oneColumn = new ArrayList<String>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			if (!line.isEmpty()) {

				String[] chunks = line.trim().split(delimiter);
				oneColumn.add(chunks[colNumber]);
			}
		}
		sc.close();
		return oneColumn;
	}
	
	
	public ArrayList<String> readColX_String_SemiColon(String fileName, int colNumber) {

		FileReader fin = null;
		try {
			fin = new FileReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);
		sc.nextLine(); // because we have headers
		ArrayList<String> oneColumn = new ArrayList<String>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			if (!line.isEmpty()) {

				String[] chunks = line.trim().split(";");
//				oneColumn.add(chunks[colNumber]);
//				oneColumn.add(chunks[colNumber].replaceAll("^\"|\"$", ""));
				oneColumn.add(chunks[colNumber].replaceAll("\"",""));
			}
		}
		sc.close();
		return oneColumn;
	}
	
	public ArrayList<Integer> readColX_Int(String fileName, int colNumber) {

		FileReader fin = null;
		try {
			fin = new FileReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);
		
		ArrayList<Integer> oneColumn = new ArrayList<Integer>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			if (!line.isEmpty()) {
				String[] chunks = line.trim().split("\t");
				oneColumn.add(Integer.parseInt(chunks[colNumber]));
				
			}
		}
		sc.close();
		return oneColumn;
	}
	
	
	public ArrayList<Integer> readColX_Int_Delimiter(String fileName, int colNumber, String delimiter, boolean header) {

		FileReader fin = null;
		try {
			fin = new FileReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);
		
		if(header){
			sc.nextLine();
		}
		ArrayList<Integer> oneColumn = new ArrayList<Integer>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			if (!line.isEmpty()) {
				String[] chunks = line.split(delimiter);

//				if(! chunks[colNumber].equalsIgnoreCase("Not found in map file")){
				oneColumn.add(Integer.parseInt(chunks[colNumber]));
//				}
			}
		}
		sc.close();
		return oneColumn;
	}
	
	public ArrayList<Double> readColX_double(String fileName, int colNumber) {

		FileReader fin = null;
		try {
			fin = new FileReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);
		ArrayList<Double> oneColumn = new ArrayList<Double>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			if (!line.isEmpty()) {
				String[] chunks = line.split("\t");
				oneColumn.add(Double.parseDouble(chunks[colNumber]));
			}
		}
		sc.close();
		return oneColumn;
	}
	
	
	
	public ArrayList<Double> readColX_double_Delimiter(String fileName, int colNumber, String delimiter, boolean header) {

		FileReader fin = null;
		try {
			fin = new FileReader(fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fin);
		ArrayList<Double> oneColumn = new ArrayList<Double>();
		
		if(header){
			sc.nextLine();
		}

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			if (!line.isEmpty()) {
				String[] chunks = line.split(delimiter);

				oneColumn.add(Double.parseDouble(chunks[colNumber]));
			}
		}
		sc.close();
		return oneColumn;
	}
	
	
	
	/** Search Functions**/
	
	public boolean searchIntList_boolean(ArrayList<Integer> arl_int, int prob){
		
		boolean result = false;
		
		for(int n: arl_int){
			if(n == prob){ result = true;}
		}
		
		return result;
	}
	
	
	
	public boolean searchIntArray_boolean(int[] a, int prob) {

		boolean b = false;

		for (int i : a) {
			if (i == prob) {
				b = true;
			}
		}

		return b;
	}
	
	public boolean searchDoubleArray_boolean(double[] a, double prob) {

		for (double d : a) {
			
			if (Math.abs(d - prob) < 0.000000001) {
			  return true;
			}
		}

		return false;
	}
	
	public boolean searchListInteger_boolean(int prob, ArrayList<Integer> myList) {

		for (int i : myList) {
			if (i == prob) {
				return true;
			}
		}

		return false;
	}
	
	public boolean searchListString_boolean(String prob, ArrayList<String> myList) {

		boolean b = false;

		for (String i : myList) {

			if (i.equals(prob)) {
				b = true;
			}
		}
		return b;

	}
	
	public boolean searchListString_boolean_ignorCase(String prob, ArrayList<String> myList) {

		boolean b = false;

		for (String i : myList) {

			if (i.equalsIgnoreCase(prob)) {
				b = true;
			}
		}
		return b;

	}
	
	public int searchListString_index(String prob, List<String> myList){
		return myList.indexOf(prob);
	}

	
	public double calcMeanArray_Int(int [] myArray){
		
		double mean  =0;
		
		for(int i: myArray){
			
			mean += i;
		}
		
		return (mean/myArray.length);
		
	}
	
	public double[] convertIntToDoubleArray(int [] a) {
		double[] b = new double[a.length];

		for (int p = 0; p < a.length; p++) {
			b[p] = a[p] * 1.0;
		}
		
		return b;
	}

}
