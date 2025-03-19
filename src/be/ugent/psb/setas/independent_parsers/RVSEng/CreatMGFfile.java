package be.ugent.psb.setas.independent_parsers.RVSEng;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class CreatMGFfile {

	/* */
	public List<List<Integer>> read_all(String filename) {

		FileReader fin = null;
		try {
			fin = new FileReader(filename);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		/* The first line is a header: */
		Scanner sc = new Scanner(fin);
		sc.nextLine();

		List<List<Integer>> gfCounts = new ArrayList<List<Integer>>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();
			String[] chunks = line.split("\t");

			List<Integer> geneCounts = new ArrayList<Integer>();

			// column 0: GF-ID , 1: num of genes, 2: num of Sp
			for (int i = 3; i < chunks.length; i++) {
				int parsed = Integer.parseInt(chunks[i]);
				geneCounts.add(parsed);
			}
			gfCounts.add(geneCounts);
		}
		sc.close();
		return gfCounts;
	}
	
	public static void main(String[] args) {

		CreatMGFfile crtMGF = new CreatMGFfile();
		CommonFunctions cmf = new CommonFunctions();
		
//		String path1 = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/combinedOutput/comSortLam_37speMGCF5_CPMpval_9178coreGF";
		String path1="/home/setas/Desktop/setas/Project1/Results/CompareRankings/newCombinedOutput/combinedOutputFiles_rawSorted/SortedLambda/37spe/comOut_37spe_Tau_Mon2_Musa3_SortedLambda";
//		String path1 = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-28Eudicots/combinedOutput/comOut_28Eud_Prun37_9178OrthGF_SortLam";
		List<List<String>> comOutput = cmf.readMapFile(path1);
		List<String> gfIDsInCombOutFile = cmf.readColX_String(path1, 0);
		
		String path2 = "/home/setas/git/StochasticBD/src/Files/GeneFamilyCounts/37spe-MGCF5-9178core-OrderNewickTree";	
		List<String> gfIDs_GFcountFile = cmf.readColX_String(path2, 0); //has a header
		List<List<Integer>> gfCounts = crtMGF.read_all(path2);
		
		
//		String path3= "/home/setas/Desktop/moreThan10WGDsRvsEng";
//		List<String> gfIDsMoreThan11wgd = cmf.read1ColFile_String(path3) ;
		
//		int numGfsMoreThan11wgd = gfIDsMoreThan11wgd.size();
		
//		for(int i=8486; i< 9178; i++){ //last 692 GFs
//		for(int i=1; i<= 970; i++){		
//		for(int i= 709; i< 8486; i++){
		
//		for(int i=0; i<numGfsMoreThan11wgd; i++){
		for(int i=0; i<9178; i++){
			
//			String gfIDProb = gfIDsMoreThan11wgd.get(i);
			String gfIDProb = gfIDs_GFcountFile.get(i+1);
			int indexInCombOutput = cmf.searchListString_index(gfIDProb, gfIDsInCombOutFile);
			
			List<String> currentLine = comOutput.get(indexInCombOutput);
			String gfID = currentLine.get(0);
			double rStar = Double.parseDouble(currentLine.get(1));
			double lambdaStar = Double.parseDouble(currentLine.get(2));
			double loglkStar = Double.parseDouble(currentLine.get(3));
//			double pvStar = Double.parseDouble(currentLine.get(4));
			
//			System.out.print(gfID+"\t"+rStar+"\t"+lambdaStar+"\t"+loglkStar+"\t"+pvStar);
			System.out.print(gfID+"\t"+rStar+"\t"+lambdaStar+"\t"+loglkStar);
			
			int index = cmf.searchListString_index(gfID, gfIDs_GFcountFile); 
			
			for(int j=0; j< 37 ; j++){
//			for(int j=19; j< 27 ; j++){
				
			System.out.print("\t"+gfCounts.get(index-1).get(j)); //the indexes would be shifted by one due to header
			
			}
			
			System.out.println();
			
			
		}
		
	}
	
	
	
}
