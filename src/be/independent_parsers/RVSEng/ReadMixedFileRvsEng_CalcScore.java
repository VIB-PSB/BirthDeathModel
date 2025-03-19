package be.ugent.psb.setas.independent_parsers.RVSEng;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class ReadMixedFileRvsEng_CalcScore {

	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();

		String path1 = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/Tair/GFid_AthaGene_InOrder";
		ArrayList<String> gfIDs = cmf.readColX_String(path1, 0);
		ArrayList<String> AthaGenes = cmf.readColX_String(path1, 1);

		String path = "/home/setas/Desktop/summaryRvsEng_AthaGenes";

		List<List<String>> map = cmf.readMapFile(path);
		// String myRegex = "ORTHO";

		// List<List<String>> higherLogLowerLambda = new
		// ArrayList<List<String>>();

		DecimalFormat df = new DecimalFormat("0.0000");

		for (int i =1; i < map.size(); i++) { // To Skip the header start at
												// i=1;
			int score = 0;
			double sumDiffLoglk = 0;
			List<String> currentGFrecord = map.get(i);

			if (!currentGFrecord.get(4).isEmpty()) {

//				double refLoglk = Double.parseDouble(currentGFrecord.get(4));
				double refLoglk = Double.parseDouble(currentGFrecord.get(12));

//				for (int j = 5; j <= 22; j++) { // 18 WGDs
				for (int j = 13; j <= 30; j++) { // 18 WGDs

					if (!currentGFrecord.get(j).isEmpty()) {
						double loglkRm1WGD = Double.parseDouble(currentGFrecord
								.get(j));

//						sumDiffLoglk += (loglkRm1WGD - refLoglk);

						if (loglkRm1WGD < refLoglk) {
							score += 1;
						}
					}
				}
			}

//			System.out.println(sumDiffLoglk*1.00/18);
			System.out.println(score);

		}

//		 for (int i = 0; i < map.size(); i++) { // To Skip the header, i=1
//		
//		 List<String> currentGFrecord = map.get(i);
//		 int rank = Integer.parseInt(currentGFrecord.get(0));
//		 String gfIDProb = currentGFrecord.get(1);
//		
//		 int maxNumAthaInMGFs = 13;
//		 String[] AthaGenesArray = new String[maxNumAthaInMGFs];
//		
//		 System.out.print(rank+"\t"+gfIDProb + "\t");
//		
//		 int counter = 0;
//		
//		 for (int lineNum = 0; lineNum < 13042; lineNum++) {
//		
//		 if (gfIDs.get(lineNum).equals(gfIDProb)) {
//		
//		 AthaGenesArray[counter] = AthaGenes.get(lineNum);
//		 counter++;
//		 }
//		 }
//		
//		 for (int l = 0; l < AthaGenesArray.length; l++) {
//		
//		 System.out.print(AthaGenesArray[l] + "\t");
//		 }
//		 for (int j = 2; j < currentGFrecord.size(); j++) { // indexes 0 and 1 are already printed
//		
//		 if(!currentGFrecord.isEmpty()){
//		 System.out.print(currentGFrecord.get(j) + "\t");}
//		
//		 else{
//		 System.out.print(""+"\t");
//		 }
//		 }
//		
//		 System.out.println();
//		
//		 }
	}
}
