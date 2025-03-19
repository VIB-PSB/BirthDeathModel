package be.ugent.psb.setas.independent_parsers;

import java.util.ArrayList;
import java.util.List;

public class GFsizeDistribution {

	 public int sumIntsInList(List<Integer> lsint){
	 int sum =0;
	 for(int i: lsint){sum +=i; }
	 return sum;}

	public ArrayList<Integer> getDiffTwoIntList(List<Integer> lsa,
			ArrayList<Integer> lsb) {

		ArrayList<Integer> diff = new ArrayList<Integer>(lsa.size());

		if (lsa.size() != lsb.size()) {
			System.out.println("error: size mismatch");}
		for (int i = 0; i < lsa.size(); i++) {
			int df = lsa.get(i) - lsb.get(i);
			diff.add(i, df);}
		return diff;
	}
	
	public ArrayList<Integer> getAbsDiffTwoIntList(ArrayList<Integer> lsa,
			ArrayList<Integer> lsb) {

		ArrayList<Integer> diff = new ArrayList<Integer>(lsa.size());

		if (lsa.size() != lsb.size()) {
			System.out.println("error: size mismatch");}
		for (int i = 0; i < lsa.size(); i++) {
			int df = Math.abs(lsa.get(i) - lsb.get(i));
			diff.add(i, df);}
		return diff;
	}

	public ArrayList<Integer> multiplyByXList(ArrayList<Integer> lsa, int x) {
		ArrayList<Integer> doubled = new ArrayList<Integer>(lsa.size());

		for (int i = 0; i < lsa.size(); i++) {

			int doubledValue = x * lsa.get(i);
			doubled.add(i, doubledValue);
		}
		return doubled;
	}

	public int maxDiffList(ArrayList<Integer> diff_lsa) {

		int max = diff_lsa.get(0);
		for (int i = 1; i < diff_lsa.size(); i++) {
			if (Math.abs(diff_lsa.get(i)) > Math.abs(max)) {
				max = diff_lsa.get(i);
			}
		}
		return max; // in case -7 and 7, we return one of the two, doesn't matter
	}

	public ArrayList<Integer> returnMaxIndexes(ArrayList<Integer> lsa, int max) {

		ArrayList<Integer> maxIndexes = new ArrayList<Integer>();
		int maximumValue = maxDiffList(lsa);

		for (int i = 0; i < lsa.size(); i++) {
			if (Math.abs(lsa.get(i)) == Math.abs(maximumValue)) {
				maxIndexes.add(i);
			}
		}
		return maxIndexes; // here we return both species with 7 and -7 indexes,
							// then we use these indexes to return both true
							// maximum values
	}

//	public static void main(String[] args) {
//
//		CommonFunctions cmmFunct = new CommonFunctions();
//		GFsizeDistribution gfsizedist = new GFsizeDistribution();
//
//		 String path1 ="/home/setas/git/IndependentParsers/Independent Parsers/src/files/37AngioSperms/baseFiles/9178coreGF_inOrderLam";
//		 List<String> GF_ids_order = cmmFunct.read1ColFile_String(path1);
//
////		 String path2 = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/combinedOutput/comSortLam_37speMGCF5_CPMpval_9178coreGF";
//		String path2 = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-28Eudicots/combinedOutput/comOut_28Eud_Prun37_9178OrthGF_SortLam";
////		 String path2="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-8Monocots/9178_8Monocots_sorted_lambda";
//		
//		ArrayList<String> gfIDs_combinedOutput = cmmFunct.readColX_String(path2, 0);
//		ArrayList<Integer> optimalRootSizes = cmmFunct.readColX_Int(path2, 1);
//
//		/* all these lists are prepared in the same order as gfcounts: newick tree (post)order */
////		String idealProfilePath = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ComparisonsWithIdealprofiles/IdealProfiles_R1";
//		String idealProfilePath = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-28Eudicots/ComparisonsWithIdealprofiles/idealProfiles_Eud.txt";
////		String idealProfilePath="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-8Monocots/ComparisonsWithIdealprofiles/idealProfiles_Monocots";
//
//		ArrayList<String> speInOrder = cmmFunct.readColX_String(
//				idealProfilePath, 0);
//		ArrayList<Integer> idealProfileR1 = cmmFunct.readColX_Int(
//				idealProfilePath, 1);
//		// ArrayList<Integer> idealProfileR1_Musa3 =
//		// cmmFunct.readColX_Int(idealprofPath, 2);
//		// ArrayList<Integer> idealProfileR1_MonocotTau =
//		// cmmFunct.readColX_Int(idealprofPath, 3);
//		// ArrayList<Integer> idealProfileR1_MonocotTau_Musa3 =
//		// cmmFunct.readColX_Int(idealprofPath, 4);
//		// ArrayList<Integer> idealProfileR1_AllAngiosperm =
//		// cmmFunct.readColX_Int(idealprofPath, 5);
//		// ArrayList<Integer> idealProfileR1_AllAngiosperm_Musa3 =
//		// cmmFunct.readColX_Int(idealprofPath, 6);
//		// ArrayList<Integer> idealProfileR1_AllAngiosperm_MonocotTau =
//		// cmmFunct.readColX_Int(idealprofPath, 7);
//		// ArrayList<Integer> idealProfileR1_AllAngiosperm_Monocottau_Musa3 =
//		// cmmFunct.readColX_Int(idealprofPath, 8);
//
//		ReadGFcountsFile rgf = new ReadGFcountsFile();
//		List<List<Integer>> gfCounts = rgf.read_all("/home/setas/git/StochasticBD/src/Files/GeneFamilyCounts/37spe-MGCF5-9178core-OrderNewickTree");
//		ArrayList<String> gfIDs_geneCountFile = rgf.getGfIDs();
////		DecimalFormat df = new DecimalFormat("0.000");
//
////		for (int i = 0; i < GF_ids_order.size(); i++) { // to be able to have all output files in one order
////
////			String gfID_current = GF_ids_order.get(i);
////			
////			int indexInComOutFile = cmmFunct.searchListString_index(gfID_current, gfIDs_combinedOutput);
////			int rootSize_current = optimalRootSizes.get(indexInComOutFile);
//////			double lambda_current = optimalLambdas.get(i);
////
////			for (int j = 0; j < gfIDs_geneCountFile.size(); j++) {
////
////				if (gfIDs_geneCountFile.get(j).equalsIgnoreCase(gfID_current)) {
////
////					// System.out.print(gfID_current + "\t" + rootSize_current
////					// + "\t" + df.format(lambda_current) );
////
////					System.out.print(gfID_current);
////					List<Integer> geneCountProfile_current = gfCounts.get(j);
//////					List<Integer> geneCountProfile_current_2 = new ArrayList<Integer>(28);
//////					
//////					for(int m=0; m<28; m++){ //for Eudicots
//////						geneCountProfile_current_2.add(m,geneCountProfile_current.get(m));	
//////					}
////					
//////					for(int m=28;m<36;m++){//For monocots
//////						geneCountProfile_current_2.add(m-28,geneCountProfile_current.get(m));
//////					}
////					//
////					// for(int count : geneCountProfile_current){
////					// System.out.print("\t"+ count);
////					// }
////					ArrayList<Integer> difference_idealOptimalR = gfsizedist
////							.getDiffTwoIntList(geneCountProfile_current,
////									gfsizedist.multiplyByXList(idealProfileR1,
////											rootSize_current));
////
////					int maxValueDiff = gfsizedist
////							.maxDiffList(difference_idealOptimalR);
////					// System.out.print(maxValueDiff + "\t"); // this does not
////					// distinguish between -x and x, so we use indexes
////
////					List<Integer> indexesOfMaxDiff = gfsizedist
////							.returnMaxIndexes(difference_idealOptimalR,
////									maxValueDiff);
////
////					for (int foundIndex : indexesOfMaxDiff) {
////						
////						System.out.print("\t" + speInOrder.get(foundIndex)
////								+ "\t"
////								+ difference_idealOptimalR.get(foundIndex));
////						
//////						System.out.print("\t" + speInOrder.get(foundIndex)
//////								+ "\t"
//////								+ difference_idealOptimalR.get(foundIndex));
////					}
////
//////					for (int k = 0; k < difference_idealOptimalR.size(); k++) {
//////						System.out
//////								.print("\t" + difference_idealOptimalR.get(k));
//////					}
////				}
////			}
////			System.out.println();
////
////		}
//
//		 int [] numOfGenesPerFamily = new int[GF_ids_order.size()];
//		 String path="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/AllGenesAllspeAllGFs30PlazaIDs";
//		 ArrayList<String> myGfIDs = cmmFunct.readColX_String(path,0);
//		 ArrayList<String> myGeneIDs = cmmFunct.readColX_String(path,1);
//		
//		 String speciesTag="BR";
//		 for(int i=0; i<GF_ids_order.size(); i++){
//		
//		 String gfID_current = GF_ids_order.get(i);
//		 int numOfGenes_temp= 0;
//		
//		 for(int j=0; j<myGfIDs.size();j++){
//		
//		 if(myGfIDs.get(j).equalsIgnoreCase(gfID_current) &&
//		 myGeneIDs.get(j).split(speciesTag).length>1 ){
//		 // && myGeneIDs.get(j).split("ATR").length <= 1
//		
//		 numOfGenes_temp+=1;
//		 }
//		 }
//		
//		 numOfGenesPerFamily[i] = numOfGenes_temp;
//		
//		 System.out.println(gfID_current+"\t"+numOfGenes_temp);
//		
//		 }
//
//	}
}
