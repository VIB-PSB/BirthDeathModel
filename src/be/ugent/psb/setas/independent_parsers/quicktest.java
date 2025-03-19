package be.ugent.psb.setas.independent_parsers;

import java.util.ArrayList;

public class quicktest {
	

	public static String findRep(ArrayList<String> arln) {
		String d = "";
		for (int i = 1; i < arln.size(); i++) {

			String s1= arln.get(i);
			String s2= arln.get(i-1);
			
			double ps1= Double.parseDouble(s1);
			double ps2 = Double.parseDouble(s2);
			
			if (Math.abs(ps1-ps2) < 1.0) {

				System.out.println(ps1);
		
			}

		}

		return d;
	}
	
	
	public static boolean ifvalueIsInArrayList(ArrayList<String> myar, double prob) {
		
		boolean b = false;
		
		for(String s: myar) {
			
			double ps = Double.parseDouble(s);
			
			if(Math.abs(ps-prob)<1.0) {
				b=true;
			}
		}
		
		return b;
	}

	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();

//		String path1 = "/home/setas/Desktop/lambdaRanking_BED_negatives";
//		String path2 = "/home/setas/Desktop/lambdaRanking_BED_positives";
		
		String path1 ="/home/setas/Desktop/lambdaRanking_UEDrepressed_positive";
		String path2="/home/setas/Desktop/lambdaRanking_UEDrepressed_negative";

		ArrayList<String> pos = cmf.read1ColFile_String(path2);
		ArrayList<String> neg = cmf.read1ColFile_String(path1);

//		 for(String p: pos) {
//		
//		 double positive = Double.parseDouble(p);
//		
//		 for(String n: neg) {
//		
//		 double negative = Double.parseDouble(n);
//		
//		 if(Math.abs(positive-negative) < 1.0) {
//		
//		 System.out.println(p+"\t"+n);
//		 }
//		 }
//		
//		 }
		
//		System.out.println(findRep(pos));
		
		double[] allRanks = new double [9178];
		
		for(int i=0; i< 9178;i++) {
			
			allRanks[i]=(i+1)*1.0;
		}
		
		for(double d: allRanks) {
			
			if(!ifvalueIsInArrayList(pos,d) && !ifvalueIsInArrayList(neg,d)) {
				
				System.out.println(d);
				
			}
		}
	

	}

}
