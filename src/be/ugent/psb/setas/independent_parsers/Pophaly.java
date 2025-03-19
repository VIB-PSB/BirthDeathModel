package be.ugent.psb.setas.independent_parsers;

import java.util.ArrayList;

public class Pophaly {

	public int[] creatComplemantoryRanks(int[] positiveRanks, int numOfGeneFamilies) {

		CommonFunctions cmf = new CommonFunctions();
		int[] b = new int[numOfGeneFamilies - (positiveRanks.length)+10];
		int index = 0;
		for (int rankValue = 1; rankValue <= numOfGeneFamilies; rankValue++) {
			if (!cmf.searchIntArray_boolean(positiveRanks, rankValue)) {
				b[index] = rankValue;
				index = index + 1;
			}
		}
		return b;
	}

	public int[] convert_IntList_to_intarray(ArrayList<Integer> ari) {

		int[] converted = new int[ari.size()];

		for (int i = 0; i < ari.size(); i++) {

			converted[i] = ari.get(i).intValue();
		}
		return converted;
	}

	public double[] convert_int_double_array(int[] a) {

		double[] b = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			b[i] = a[i] * 1.0;
		}
		return b;
	}

//	public static void main(String[] args) {
//
//		CommonFunctions cmf = new CommonFunctions();
//
//		/** To creat GF Gene list of Zmay genes for comparing results with Pophaly */
//		// String path_GF_Genes
//		// ="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/AllGenes_AllGFs_37pePlazaIDs";
//		//
//		// List<List<String>> map = cmf.readMapFile(path_GF_Genes);
//		// for(int i=0; i< map.size();i++){
//		//
//		// List<String> line = map.get(i);
//		//
//		// if(line.get(1).split("ZM").length >1 ){
//		//
//		// System.out.println(line.get(0)+"\t"+line.get(1));
//		// }
//		// }
//
//		// String path_GF_ZM="/home/setas/Desktop/GF_ZMgenes";
//		// String path_MibelFile="/home/setas/Desktop/maize_id_mapping.txt";
//		//
//		// List<List<String>> map_GF_ZM = cmf.readMapFile(path_GF_ZM);
//		// List<List<String>> map_mibel = cmf.readMapFile(path_MibelFile);
//		//
//		// for(int i=0;i<map_GF_ZM.size();i++){
//		//
//		// List<String> line = map_GF_ZM.get(i);
//		//
//		// for(int j=0; j<map_mibel.size();j++){
//		//
//		// List<String> newLine = map_mibel.get(j);
//		//
//		// if(line.get(1).equalsIgnoreCase(newLine.get(0))){
//		//
//		// System.out.println(line.get(0)+"\t"+line.get(1)+"\t"+newLine.get(1)+"\t"+newLine.get(2));
//		// }
//		// }
//		// }
//
//		/** end of making input file. Feed this to the next part **/
//
//		/** Add a column BED or UED to the previous file **/
//		// String path_GF_ZM_GMZ =
//		// "/home/setas/Desktop/setas/setas_stmae/Reviews/Pophaly/GF_ZMgenes_GMZids";
//		////// String Pophaly
//		// ="/home/setas/Desktop/setas/setas_stmae/Reviews/Pophaly/subfunctionalized_BEDgenes";
//		// String
//		// Pophaly="/home/setas/Desktop/setas/setas_stmae/Reviews/Pophaly/input/dominant_UED";
//		////
//		// List<List<String>> map_GF_ZM = cmf.readMapFile(path_GF_ZM_GMZ);
//		// List<List<String>> map_pophaly = cmf.readMapFile(Pophaly);
//		//
//		// for(int i=0; i<map_GF_ZM.size();i++){
//		//
//		// List<String> line = map_GF_ZM.get(i);
//		// boolean b = false;
//		//
//		// for(int j=0;j< map_pophaly.size();j++){
//		//
//		// List<String> pophalyLine = map_pophaly.get(j);
//		//
//		// if(line.get(2).equalsIgnoreCase(pophalyLine.get(0))||
//		// line.get(3).equalsIgnoreCase(pophalyLine.get(0))){
//		//
//		// b = true;
//		// System.out.println(line.get(0)+"\t"+line.get(1)+"\t"+line.get(2)+"\t"+pophalyLine.get(0));
//		// }
//		//
//		// }
//		//
//		//// if(b==false){ //to print all BED and non BED
//		//// System.out.println(line.get(0)+"\t"+line.get(1)+"\t"+line.get(2));
//		//// }
//		//
//		// }
//		/** end of making the final input file for comparing to Pophaly **/
//
//		Pophaly pof = new Pophaly();
//		String pathcombinedOutputFile = "/home/setas/Desktop/setas/setas_stmae/Reviews/CombindeOutputCompletion/TableS1_PlantCell_noHeader.txt";
//		String path_BEDGenes = "/home/setas/Desktop/setas/setas_stmae/Reviews/Pophaly/GF_ZMgenes_GMZids_BEDgenes";
//		String path_dominantUED = "/home/setas/Desktop/setas/setas_stmae/Reviews/Pophaly/GF_ZMgenes_GMZids_DominantUED";
//		String path_repressedUED = "/home/setas/Desktop/setas/setas_stmae/Reviews/Pophaly/GF_ZMgenes_GMZids_UEDrepressed";
//
//		List<List<String>> map_combinedOutput = cmf.readMapFile(pathcombinedOutputFile);
//		List<String> gfIDs_BED = cmf.readColX_String(path_BEDGenes, 0);
//		List<String> gfIDs_Dominant_UED = cmf.readColX_String(path_dominantUED, 0);
//		List<String> gfIDs_Rep_UED = cmf.readColX_String(path_repressedUED, 0);
//
//		 ArrayList<Double> positiveRanks = new ArrayList<Double>();
////		ArrayList<Integer> posRanks = new ArrayList<Integer>();
//		CreatRankVectorWilcoxonTest crtRank = new CreatRankVectorWilcoxonTest(cmf);
//
//		for (int j = 0; j < map_combinedOutput.size(); j++) {
//
//			List<String> combLine = map_combinedOutput.get(j);
//			String current_gfID = combLine.get(1);
//
//			String category = "NA";
//			 if (gfIDs_BED.contains(current_gfID)) {
//			 category = "BED";
//			 positiveRanks.add(Double.parseDouble(combLine.get(7))); // lambdaRank
////			 posRanks.add(Integer.parseInt(combLine.get(0))); //rankInteger_first column
//			 }
//
////			if (gfIDs_Dominant_UED.contains(current_gfID)) {
////				category = "UED_dominant";
////				positiveRanks.add(Double.parseDouble(combLine.get(7))); // lambdaRank
//////				positiveRanks.add(Integer.parseInt(combLine.get(0)));
////			}
//			
////			 if(gfIDs_Rep_UED.contains(current_gfID)){
////			  category="UED_repressed";
////			  
//////			  if(!positiveRanks.contains(Double.parseDouble(combLine.get(7)))) {
////		      positiveRanks.add(Double.parseDouble(combLine.get(7))); // lambdaRank
//////			  positiveRanks.add(Integer.parseInt(combLine.get(0)));
//////			  }
////			 }
//		}	
//		ArrayList<Integer> positiveRanks2 = new ArrayList<Integer>();
//		for (double r : positiveRanks) {
////			System.out.println((int) (Math.round(r)));
//			positiveRanks2.add((int) (Math.round(r)));
//		}
//		// When working with integer ranks
////		int[] positives = pof.convert_IntList_to_intarray(posRanks);
//		int[] positives = pof.convert_IntList_to_intarray(positiveRanks2);
//		int[] negatives = pof.creatComplemantoryRanks(positives, 9178);
//
//		double[] pos = pof.convert_int_double_array(positives);
//		double[] neg = pof.convert_int_double_array(negatives);
//
////		System.out.println("UED repressed");
////		for(double p: pos){System.out.println(p);}
////		for (double n : neg) {System.out.println(n);}
//
//		 if (pos.length >= 2 && neg.length >= 2) {
//		 jsc.independentsamples.MannWhitneyTest man_less = new MannWhitneyTest(pos,
//		 neg, H1.NOT_EQUAL);
//		 double U = man_less.getTestStatistic();
//		 double z = man_less.getZ();
//		 double pval = man_less.approxSP();
//		 double median_pos = CreatRankVectorWilcoxonTest.calcMedian(pos);
//		 double median_neg = CreatRankVectorWilcoxonTest.calcMedian(neg);
//	
//		 System.out.println(U+"\t"+z+"\t"+pval+"\t"+median_pos+"\t"+median_neg);}
//
//		/*** OLD MWU ***/
//		// MannWhitneyUTest test = new MannWhitneyUTest();
//		// double[] us = crtRank.calculateUxUy(positives, negatives);
//		// System.out.print(us[0] + "\t" + us[1] + "\t"); // prints Ux and Uy
//		// double pvalue = test.mannWhitneyUTest(pos, neg);
//		// double u = test.mannWhitneyU(pos, neg);
//		// System.out.print(u + "\t" + pvalue + "\n");
//	}
}
