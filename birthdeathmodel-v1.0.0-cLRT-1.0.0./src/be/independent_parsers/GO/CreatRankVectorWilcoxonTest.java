package be.ugent.psb.setas.independent_parsers.GO;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;
import be.ugent.psb.setas.independent_parsers.Pophaly;
import jsc.independentsamples.MannWhitneyMedianDifferenceCI;
import jsc.independentsamples.MannWhitneyTest;
import jsc.tests.H1;

public class CreatRankVectorWilcoxonTest {

	private CommonFunctions cmmFunct;

	public CreatRankVectorWilcoxonTest(CommonFunctions cmf) {
		this.cmmFunct = cmf;

	}

	public ArrayList<String> read_String_Column_EqualSignDelimiter(String fileName, int colNumber) {

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
				String[] chunks = line.split(" = ");
				oneColumn.add(chunks[colNumber]);
			}
		}
		sc.close();
		return oneColumn;
	}

	public int[] search_GoID(int GoID_prob, ArrayList<Integer> GoIDs, ArrayList<Integer> ranks, int lastRank,
			ArrayList<String> GFids) {

		int[] found = new int[lastRank];

		for (int i = 0; i < GoIDs.size(); i++) {

			if (GoIDs.get(i) == GoID_prob) {

				int rank = ranks.get(i);
				/**
				 * we don't have rank 0, so index=0 corresponds to no family and will be always
				 * 0 by default which can result in wrong data for simulation, so it is skipped
				 */
				found[rank - 1] = 1;

			}
		}
		return found;
	}

	public int[] search_GoID_returnRealRank(int GoID_prob, ArrayList<Integer> GoIDs, ArrayList<Integer> ranks) {

		ArrayList<Integer> found = new ArrayList<Integer>();

		for (int i = 0; i < GoIDs.size(); i++) {

			if (GoIDs.get(i) == GoID_prob) {

				int rank = ranks.get(i);
				found.add(rank);

			}
		}

		return cmmFunct.convertListToArray(found);
	}

	public ArrayList<Double> search_GoID_returnRealRank_double(int GoID_prob, ArrayList<Integer> GoIDs,
			ArrayList<Double> ranks) {

		ArrayList<Double> foundRanks = new ArrayList<Double>();

		for (int i = 0; i < GoIDs.size(); i++) {

			if (GoIDs.get(i) == GoID_prob) {

				double rank = ranks.get(i);
				foundRanks.add(rank);

			}
		}

		// return cmmFunct.convertListToArray_double(found);
		return foundRanks;
	}

	public ArrayList<Double> search_families_returnRank_double(String familyName, ArrayList<String> families,
			ArrayList<String> GOdesc, ArrayList<Double> ranks) {

		ArrayList<Double> foundRanks = new ArrayList<Double>();

		for (int i = 0; i < families.size(); i++) {

			// if (families.get(i).equalsIgnoreCase(familyName) ||
			// GOdesc.get(i).contains(familyName) ) { //contain matches the exact pattern
			// if(families.get(i).equalsIgnoreCase(familyName) ||
			// GOdesc.get(i).split(familyName).length > 1){ //split is case insensitive
			// if (families.get(i).equalsIgnoreCase(familyName) ||
			// Pattern.compile(Pattern.quote(familyName),
			// Pattern.CASE_INSENSITIVE).matcher(GOdesc.get(i)).find() ) {
			if (families.get(i).contains(familyName) || GOdesc.get(i).contains(familyName)) {
				double rank = ranks.get(i);
				foundRanks.add(rank);

			}
		}

		// return cmmFunct.convertListToArray_double(found);
		return foundRanks;
	}

	public ArrayList<Integer> search_families_returnRank_Int(String familyName, ArrayList<String> families,
			ArrayList<String> GOdesc, ArrayList<Integer> ranks) {

		ArrayList<Integer> foundRanks = new ArrayList<Integer>();

		for (int i = 0; i < families.size(); i++) {

			// if (families.get(i).equalsIgnoreCase(familyName) ||
			// GOdesc.get(i).contains(familyName) ) {
			// if(families.get(i).equalsIgnoreCase(familyName) ||
			// GOdesc.get(i).split(familyName).length > 1){
			// if (families.get(i).equalsIgnoreCase(familyName) ||
			// Pattern.compile(Pattern.quote(familyName),
			// Pattern.CASE_INSENSITIVE).matcher(GOdesc.get(i)).find() ) {
			if (families.get(i).contains(familyName) || GOdesc.get(i).contains(familyName)) { // contain is
																								// case-sensitive
				int rank = ranks.get(i).intValue();
				foundRanks.add(rank);

			}
		}

		// return cmmFunct.convertListToArray_double(found);
		return foundRanks;
	}

	public int[] creatComplemantoryRanks(int[] a, int numOfGeneFamilies) {

		int[] b = new int[numOfGeneFamilies - (a.length)];
		int index = 0;

		for (int i = 1; i <= numOfGeneFamilies; i++) {

			if (!cmmFunct.searchIntArray_boolean(a, i)) {

				b[index] = i;
				index++;
			}
		}

		return b;
	}

	public ArrayList<Double> creatComplemantoryRanks_double(ArrayList<Double> foundPositive, int numOfGeneFamilies) {

		int lenght = numOfGeneFamilies - foundPositive.size();
		ArrayList<Double> negativeRanks = new ArrayList<Double>(lenght);

		for (int i = 1; i <= numOfGeneFamilies; i++) {

			double rankProb = i * 1.0;
			boolean b = false;

			// if (!foundPositive.contains(rankProb)) { //if the floor or ceiling is exactly
			// the same, this doesn't work
			for (Double d : foundPositive) {

				if (Math.abs(d.doubleValue() - rankProb) < 1.0) {
					// negativeRanks.add(rankProb);
					b = true;
				}
			}

			if (b == false) {
				negativeRanks.add(rankProb);
			}
		}

		return negativeRanks;
	}

	public ArrayList<Boolean> checkNonRepetetiveGoIds(ArrayList<Integer> GoIDs) {

		ArrayList<Integer> nonRep = new ArrayList<Integer>();
		ArrayList<Boolean> nonRep_boolean = new ArrayList<Boolean>();

		nonRep.add(GoIDs.get(0));
		nonRep_boolean.add(true);

		for (int i = 1; i < GoIDs.size(); i++) {

			if (!cmmFunct.searchListInteger_boolean(GoIDs.get(i), nonRep)) {

				nonRep.add(GoIDs.get(i));
				nonRep_boolean.add(true);
			} else {
				nonRep_boolean.add(false);
			}
		}
		return nonRep_boolean;
	}

	public static double[] convertIntToDouble(int[] hasGo) {

		double[] ranks_double = new double[hasGo.length];
		int i = 0;
		for (int value : hasGo) {
			ranks_double[i] = (double) value;
			i++;
		}

		return ranks_double;
	}

	public double[] calculateUxUy(int[] x, int[] y) {

		double sumUx = 0;
		double sumUy = 0;

		for (int xi : x) {

			for (int yi : y) {

				if (xi > yi) {

					sumUx += 1;
				}

				if (yi > xi) {
					sumUy += 1;
				}

				if (xi == yi) {

					sumUx += 0.5;
					sumUy += 0.5;
				}
			}

		}

		double[] us = new double[2];
		us[0] = sumUx;
		us[1] = sumUy;

		return us;

	}

	public double[] calculateUxUy_double(double[] x, double[] y) {

		double sumUx = 0;
		double sumUy = 0;

		for (double xi : x) {

			for (double yi : y) {

				if (xi > yi) {
					sumUx += 1;
				}

				if (yi > xi) {
					sumUy += 1;
				}

				if (xi - yi < 0.0000001) {
					sumUx += 0.5;
					sumUy += 0.5;
				}
			}
		}
		double[] us = new double[2];
		us[0] = sumUx;
		us[1] = sumUy;
		return us;
	}

	public static double calcAverage(int[] ar) {

		int sum = 0;

		for (int d : ar) {

			sum += d;
		}

		double avg = sum / (ar.length);
		return avg;
	}

	public static double calcAverage_double(double[] ar) {

		double sum = 0;

		for (double d : ar) {

			sum += d;
		}

		double avg = sum / (ar.length);
		return avg;
	}

	public ArrayList<Double> returnRefinedIntValRanks(ArrayList<Double> ranks) {

		ArrayList<Double> ranks_refined_Integers = new ArrayList<Double>();

		for (int r = 0; r < ranks.size(); r++) {

			double currentRank = ranks.get(r);
			double remainder = currentRank - Math.floor(currentRank);

			if (remainder - 0.5 <= 0.0000001) {
				ranks_refined_Integers.add(Math.floor(currentRank));
			}

			else {
				ranks_refined_Integers.add(Math.ceil(currentRank));
			}
		}
		return ranks_refined_Integers;
	}

	public static double calcMedian(double[] numArray) {

		Arrays.sort(numArray);
		double median;
		if (numArray.length % 2 == 0) {
			median = ((double) numArray[numArray.length / 2] + (double) numArray[numArray.length / 2 - 1]) / 2;
		} else {
			median = (double) numArray[numArray.length / 2];
		}

		return median;
	}

	public static void main(String[] args) {

		CommonFunctions cmmFunct = new CommonFunctions();
		CreatRankVectorWilcoxonTest crtRank = new CreatRankVectorWilcoxonTest(cmmFunct);
		int numGFsWithGOs = 9178;

		// String path = args[0]; //ranked file
		// String path
		// ="/home/setas/Desktop/setas/Project1/Results/CompareRankings/PCCranks-Matlab/NewRankedGFGODesc/ranks_GFid_GOid_Desc_37spe_TauMon2Musa3";
		// ArrayList<String> GFids = cmmFunct.readColX_String(path, 1);
		// ArrayList<Integer> GOids = cmmFunct.readColX_Int(path, 2);
		// ArrayList<String> GOdesc = cmmFunct.readColX_String(path, 3);

		/****** to take the ones with pv < 0.05 ***/
		// String
		// path2="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/MannWhitneyUtest/AngioTauMon2Musa3/GO-Desc-UxUyU-PV-AdjPV001";
		// ArrayList<Integer> goIDs_pv = cmmFunct.readColX_Int(path2, 0);
		// ArrayList<String> goDesc = cmmFunct.readColX_String(path2, 1);
		//// ArrayList<Integer> nomPosRanks = cmmFunct.readColX_Int(path2, 2);
		//// ArrayList<Integer> nomNegRanks = cmmFunct.readColX_Int(path2, 3);
		//// ArrayList<Double> muP = cmmFunct.readColX_double(path2, 4);
		//// ArrayList<Double> muN = cmmFunct.readColX_double(path2, 5);
		//
		// ArrayList<Double> uP = cmmFunct.readColX_double(path2, 2);
		// ArrayList<Double> uN = cmmFunct.readColX_double(path2, 3);
		// ArrayList<Double> U = cmmFunct.readColX_double(path2, 4);
		// ArrayList<Double> pv = cmmFunct.readColX_double(path2, 5);
		// ArrayList<Double> corrected_pv =cmmFunct.readColX_double(path2,6);

		/** To subset the ones on the top or bottom or pv less than a threshold */
		// for (int i = 0; i < goIDs_pv.size(); i++) { // all above arrayLists should
		// have the same sizes
		// if(corrected_pv.get(i) < 0.01){
		// if (uP.get(i) > uN.get(i)) { // to get the GOs tending to be on top uP <uN
		//
		// System.out.println(goIDs_pv.get(i) + "\t" + goDesc.get(i)
		// + "\t" + uP.get(i) + "\t" + uN.get(i) + "\t"+U.get(i)+"\t"
		// + pv.get(i)+"\t"+corrected_pv.get(i));
		// }
		// }
		// }

		// ArrayList<Integer> nonRep_Goids =
		// cmmFunct.nonRepetetive_IntList(ArrayList<Integer>)(GOids);
		// System.out.println("non rep size"+nonRep_Goids.size()); 8292 correct
		// in the next analysis

		////// int go_pvZero = Integer.parseInt(args[2]);
		//// int go_pvZero = 6499;
		////
		// /** To calculate pvalues and corrected ones */
		// ArrayList<Integer> visitedGoids = new ArrayList<Integer>();

		// for (int go = 0; go < GOids.size(); go++) {
		//
		//// if(GOids.get(go)== go_pvZero){ /** for one specific GO **/
		//
		// /** To avoid looking for repetetive GO IDs **/
		// if (!cmmFunct.searchListInteger_boolean(GOids.get(go), visitedGoids)) {
		//
		// visitedGoids.add(GOids.get(go));
		//
		// System.out.print(GOids.get(go) + "\t" + GOdesc.get(go));
		//
		// ArrayList<Double> posRanks =
		// crtRank.search_GoID_returnRealRank_double(GOids.get(go), GOids,
		// refinedranks);
		// double[] positive_ranks = cmmFunct.convertListToArray_double(posRanks);
		//// int lenx1 = positive_ranks.length;
		//
		// double[] negative_ranks =
		// cmmFunct.convertListToArray_double(crtRank.creatComplemantoryRanks_double(posRanks,numGFsWithGOs));
		//// int lenx2 = negative_ranks.length;
		//// int lenx = Math.max(lenx1, lenx2);
		//
		//// double mu1 =
		// CreatRankVectorWilcoxonTest.calcAverage_double(positive_ranks);
		//// double mu2 =
		// CreatRankVectorWilcoxonTest.calcAverage_double(negative_ranks);
		//// System.out.print("\t" + mu1 + "\t" + mu2);
		//
		// /*
		// * To print positive and negative rank vectors: for GO
		// * categories with pv=0 to be recalculated by Matlab
		// */
		//// for (int i = 0; i < lenx; i++) {
		////
		//// if (i >= lenx1 && i<lenx2 && negative_ranks[i]!=0) {
		////
		//// System.out.println("\t" + negative_ranks[i]);}
		////
		//// else if (i >= lenx2 && i<lenx1 && positive_ranks[i] !=0) {
		////
		//// System.out.println(positive_ranks[i]+ "\t");}
		////
		//// else if(i < lenx1 && i< lenx2){
		//// System.out.println(positive_ranks[i] + "\t"
		//// + negative_ranks[i]); }
		//// }
		////
		//// /** Mann-WhitneyU be.ugent.psb.setas.bdmodel.test */
		// MannWhitneyUTest test = new MannWhitneyUTest();
		//
		// double[] positives = convertIntToDouble(positive_ranks);
		//// double[] negatives = convertIntToDouble(negative_ranks);
		//// double[] us = crtRank.calculateUxUy(positives, negatives);
		//
		// double[] us = crtRank.calculateUxUy_double(positive_ranks, negative_ranks);
		// System.out.print("\t" + us[0] + "\t" + us[1]); //prints Ux and Uy
		//
		// double pvalue = test.mannWhitneyUTest(positive_ranks, negative_ranks);
		// double u = test.mannWhitneyU(positive_ranks, negative_ranks);
		//
		// System.out.print("\t" + u + "\t" + pvalue + "\n");
		//
		// }
		// }
		// }

		/*
		 * To find which gene family ids are not in the .anno file and why
		 * ArrayList<String> GFids_original = cmmFunct.readOneColumnFile(
		 * "/home/setas/git/IndependentParsers/Independent Parsers/src/files/GFids_notInBiNGOAnno"
		 * ); ArrayList<String> GFids_AnnoFile = crtRank.readOneColumnFile(
		 * "/home/setas/git/IndependentParsers/Independent Parsers/src/files/GFids_notInBiNGOinputAnno"
		 * ); // ArrayList<String> GFidsAnno_nonRep =
		 * crtRank.getNonRepetetiveGFIds(GFids_AnnoFile); //
		 * System.out.println(GFidsAnno_nonRep.size()); //9111
		 * 
		 * for(String gfId : GFids_original){ if(! crtRank.searchListString(gfId,
		 * GFids_AnnoFile)){ System.out.println(gfId); } }
		 */

		// /** For the reviews Freeling and Pophaly **/
		// String
		// path="/home/setas/Desktop/setas/setas_stmae/Reviews/Freeling/Freeling2008/comOutput_Freeling2008_Families.txt";
		//
		// ArrayList<String> GO_desc = cmmFunct.readColX_String(path, 2);
		// ArrayList<Integer> ranks = cmmFunct.readColX_Int(path, 0);
		// ArrayList<String> families = cmmFunct.readColX_String(path, 10);
		//
		//// ArrayList<String> visitedFamilies = new ArrayList<String>();
		//// ArrayList<Double> refinedranks = crtRank.returnRefinedIntValRanks(ranks);
		// ArrayList<Integer> refinedranks = ranks;
		//
		// for (int i = 0; i < families.size(); i++) { // edit gene family names to put
		// some in one group
		//
		// families.set(i, families.get(i).trim());
		//
		// // if(families.get(i).split("find").length > 1){ // for Calcium
		// // families.set(i, families.get(i).split("find")[0]);
		// // }
		//
		// if (families.get(i).split("find").length > 1 &&
		// families.get(i).split("LOB").length > 1) { // for LOB
		//
		// families.set(i, families.get(i).split("find")[0].trim());
		// }
		//
		// if (families.get(i).split("sample").length > 1) { //for F box
		//
		// families.set(i, families.get(i).split("sample")[0].trim());
		// }
		// }
		//
		// String familyProb= "GRAS";
		//
		// System.out.println(familyProb+"\n");
		// ArrayList<Integer> posRanks = new ArrayList<Integer>();
		//
		// Pophaly pof = new Pophaly();
		//
		//// for (int i = 1; i < families.size(); i++) {
		//
		//// if(!families.get(i).equalsIgnoreCase("NA")){
		//
		//// if (!cmmFunct.searchListString_boolean_ignorCase(families.get(i),
		// visitedFamilies)) {
		//// visitedFamilies.add(families.get(i));
		//// System.out.print(families.get(i)+"\t");
		//
		//// ArrayList<Double> posRanks =
		// crtRank.search_families_returnRank_double(families.get(i), families,
		// GO_desc,refinedranks);
		// posRanks = crtRank.search_families_returnRank_Int(familyProb, families,
		// GO_desc,refinedranks);
		//// ArrayList<Integer> posRanks =
		// crtRank.search_families_returnRank_Int(families.get(i), families,
		// GO_desc,refinedranks);
		//
		// int[] positives = pof.convert_IntList_to_intarray(posRanks);
		// int[] negatives = pof.creatComplemantoryRanks(positives, 9178);
		//
		// double[] positive_ranks = pof.convert_int_double_array(positives);
		// double[] negative_ranks = pof.convert_int_double_array(negatives);
		//
		// //// ArrayList<Double> refinedPosRanks =
		// crtRank.returnRefinedIntValRanks(posRanks);
		// //// double[] positive_ranks =
		// cmmFunct.convertListToArray_double(refinedPosRanks);
		// //// double[] negative_ranks =
		// cmmFunct.convertListToArray_double(crtRank.creatComplemantoryRanks_double(refinedPosRanks,numGFsWithGOs));
		// MannWhitneyUTest test = new MannWhitneyUTest();
		// double[] us = crtRank.calculateUxUy_double(positive_ranks, negative_ranks);
		//// System.out.print(us[0] + "\t" + us[1]+"\t"); //prints Ux and Uy
		// double pvalue = test.mannWhitneyUTest(positive_ranks, negative_ranks);
		// //wilcoxon_rankSum in R, paired=false
		// double u = test.mannWhitneyU(positive_ranks, negative_ranks);
		//// System.out.print(u + "\t" + pvalue + "\n");
		//// }
		//// }
		//// }

		// /** To repeat Go enrichment analysis **/
		// Do NOT read the ranks from these files, they are wrong!!!
//		String path = "/home/setas/Desktop/setas/setas_stmae/Reviews/GO_enrichments_review/rank_GF_GO_Desc_MusaFirst2Close_LambdaRanking_correct.txt";
		// String path
		// ="/home/setas/Desktop/setas/setas_stmae/Reviews/GO_enrichments_review/printIntRanksForCombinedRank_Musafirst2";

		String path = "C:\\Users\\setar\\git\\StochasticBD\\data\\allRanks_pv_H0Rm_csv.csv";
		ArrayList<String> GOIDs = cmmFunct.readColX_String_Delimiter(path, 3, "\t", true);
		Pophaly pof = new Pophaly();

		/** make non-repetettive GO ID list, to make the running time less **/
		ArrayList<String> nonRep_Goids = new ArrayList<String>();

		for (int i = 0; i < GOIDs.size(); i++) {

			String go_terms = GOIDs.get(i);
			String[] go_term = go_terms.split(";");

			for (String go : go_term) {
				if (!nonRep_Goids.contains(go)) {
					nonRep_Goids.add(go);
				}
			}
		}
		/***/
		ArrayList<String> gfIDs_comb = cmmFunct.readColX_String_Delimiter(path, 0, "\t", true);
		ArrayList<Double> ranks_comb2 = cmmFunct.readColX_double_Delimiter(path, 7, "\t"
				+ "", true);
		ArrayList<Integer> ranks_comb = new ArrayList<Integer>();
		for (double r : ranks_comb2) {ranks_comb.add((int) (Math.round(r)));}

		/** Calculate ranks for MWU **/

		for (int i = 0; i < nonRep_Goids.size(); i++) {

			String goID_prob = nonRep_Goids.get(i);
			ArrayList<Integer> positiveRanks = new ArrayList<Integer>();

			for (int lineInFile = 0; lineInFile < gfIDs_comb.size(); lineInFile++) {

				String[] goTerms = GOIDs.get(lineInFile).split(";");

				for (String goTerm : goTerms) {

					if (goTerm.equalsIgnoreCase(goID_prob)) {

						if (!positiveRanks.contains(ranks_comb.get(lineInFile))) {
							positiveRanks.add(ranks_comb.get(lineInFile));
						}
					}
				}
			}

			int[] positives_int = pof.convert_IntList_to_intarray(positiveRanks);
			int[] negatives_int = pof.creatComplemantoryRanks(positives_int, 9178);

			/** start of MWU **/
			double[] positives = cmmFunct.convertIntToDoubleArray(positives_int);
			double[] negatives = cmmFunct.convertIntToDoubleArray(negatives_int);

			if (positives.length >= 2 && negatives.length >= 2) {
				jsc.independentsamples.MannWhitneyTest man_less = new MannWhitneyTest(positives, negatives,
						H1.NOT_EQUAL);

//				 jsc.independentsamples.MannWhitneyMedianDifferenceCI = new
//				 MannWhitneyMedianDifferenceCI(a,b,0.95);

				double U = man_less.getTestStatistic();
				double z = man_less.getZ();
				double pval = man_less.approxSP();

				double median_pos = CreatRankVectorWilcoxonTest.calcMedian(positives);
				double median_neg = CreatRankVectorWilcoxonTest.calcMedian(negatives);

//				if (median_pos < median_neg) {
//					System.out.println(goID_prob + "\t" + median_pos + "\t" + median_neg + "\t" + U + "\t"
//							+ z + "\t" + pval + "\t" + "top");
//				}

//				if (median_pos > median_neg) {
//					System.out.println(goID_prob + "\t" + median_pos + "\t" + median_neg + "\t" + U + "\t" + z + "\t"
//							+ pval + "\t" + "bottom");
//				}
				
				System.out.println(goID_prob + "\t" + median_pos + "\t" + median_neg + "\t" + U + "\t"
						+ z + "\t" + pval + "\t" +"all");
			}

			/** end of MWU **/

			// jsc.independentsamples.MannWhitneyTest man_greater = new
			// MannWhitneyTest(a,b,H1.GREATER_THAN);
			// double U_greater = man_greater.getTestStatistic();
			// double z_greater= man_greater.getZ();
			// double pval_greater = man_greater.approxSP();
			// System.out.println(U_greater+"\t"+z_greater+"\t"+pval_greater);

			// double[] positive_ranks_double = pof.convert_int_double_array(positives);
			// double[] negative_ranks_double = pof.convert_int_double_array(negatives);
			//
			// WilcoxonRankSum wRankSum = new
			// WilcoxonRankSum(positive_ranks_double,negative_ranks_double);
			// String nullH = wRankSum.getNullHypothesis();
			//
			// System.out.println(nullH);
			// } // if goID prob
		} // for loop

	} // main
} // class
