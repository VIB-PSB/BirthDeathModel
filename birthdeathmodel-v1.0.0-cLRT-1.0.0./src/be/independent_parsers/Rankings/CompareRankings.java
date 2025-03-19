package be.ugent.psb.setas.independent_parsers.Rankings;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;



public class CompareRankings {
	
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
				
				for(int j=0;j<chunks.length;j++){
				ls.add(chunks[j]);
				}

				map.add(ls);
			}
		}
		sc.close();
		return map;
	}
	
	public List<List<String>> searchMap(List<List<String>> map, String probe) {
		
		List<List<String>> found = new ArrayList<List<String>>();
	
		for (List<String> ls : map) {
			// careful, now the probe must be in the column: [0]
			if (ls.get(0).equals(probe)) {
				
				found.add(ls);

			}
		}
		
		return found;

	}


	public static void main(String[] args) {

		CompareRankings cmpRank = new CompareRankings();
		CommonFunctions cmmFunc = new CommonFunctions();

		// ArrayList<String> gfIDs_zhenMapToSetas =
		// cmmFunc.read1ColFile_String("/home/setas/git/StochasticBD/src/Files/notNullSTmaeIDs_corrZhenOrderLam_toCompare");
		// ArrayList<String> gfIDs_eud =
		// cmmFunc.read1ColFile_String("/home/setas/git/StochasticBD/src/Files/stmaeGFIDs_corresToNotNullinZhen");
		//
		// if(gfIDs_zhenMapToSetas.size() != gfIDs_eud.size()){
		// System.out.println("Warning: The size of the two lists are not equal");
		// }
		//
		// for (int i =0; i<gfIDs_zhenMapToSetas.size();i++){
		// String gfID_zhenMapSetas = gfIDs_zhenMapToSetas.get(i);
		//
		// for(int j=0; j<gfIDs_eud.size();j++){
		//
		// String gfID_eud = gfIDs_eud.get(j);
		//
		// if(gfID_zhenMapSetas.equals(gfID_eud)){
		//
		// // if(i<= 300 && j> 300){
		// System.out.println(gfID_zhenMapSetas+"\t"+(i+1)+"\t"+(j+1));
		// // }
		// }
		// }
		// }

		// String path1 =
		// "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-8Monocots/ranks_comOut_8Monocots9178";
		String comAngPath = "/home/setas/Desktop/setas/Project1/Results/CompareRankings/ranksCorr_stmae_comOut_37Angio9178";
		String combEudPath = "/home/setas/Desktop/setas/Project1/Results/CompareRankings/ranksCorr_stmae_comOut_28Eud9178";
		String comMonPath = "/home/setas/Desktop/setas/Project1/Results/CompareRankings/ranksCorr_stmae_comOut_8Mon9178";

		// This file is sorted, so the order can be used to print ranks
		ArrayList<Integer> ranks_angio = cmmFunc.readColX_Int(comAngPath, 0);
		ArrayList<String> gfIDs_angio = cmmFunc.readColX_String(comAngPath, 1);
		ArrayList<Integer> optR_angio = cmmFunc.readColX_Int(comAngPath, 2);
		ArrayList<Double> optLambdas_angio = cmmFunc.readColX_double(comAngPath, 3);
		ArrayList<Double> opt_likelihoods_angio = cmmFunc.readColX_double(comAngPath, 4);
		ArrayList<Double> pvalues_angio = cmmFunc.readColX_double(comAngPath, 5);
		
		ArrayList<Integer> ranks_eud = cmmFunc.readColX_Int(combEudPath, 0);
		ArrayList<String> gfIDs_eud = cmmFunc.readColX_String(combEudPath, 1);
		ArrayList<Integer> optR_eud = cmmFunc.readColX_Int(combEudPath, 2);
		ArrayList<Double> optLambdas_eud = cmmFunc.readColX_double(combEudPath,3);
		ArrayList<Double> opt_likelihoods_eud = cmmFunc.readColX_double(combEudPath, 4);
		ArrayList<Double> pvalues_eud = cmmFunc.readColX_double(combEudPath, 5);


		ArrayList<Integer> ranks_monocots = cmmFunc.readColX_Int(comMonPath, 0);
		ArrayList<String> gfIDs_monocots = cmmFunc.readColX_String(comMonPath,1);
		ArrayList<Integer> optR_monocots = cmmFunc.readColX_Int(comMonPath, 2);
		ArrayList<Double> optLambdas_monocots = cmmFunc.readColX_double(comMonPath, 3);
		ArrayList<Double> opt_likelihoods_monocots = cmmFunc.readColX_double(comMonPath, 4);
		ArrayList<Double> pvalues_monocots = cmmFunc.readColX_double(comMonPath, 5);
		List<List<String>> mapGFAtha = cmpRank.readMapFile("/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/Tair/addedInforFromTair/GFid_Ath_geneDescrTair_InOrder");
		
		
//		DecimalFormat df = new DecimalFormat("0.000");

		for (int i = 0; i < 9178; i++) {

			String gfID_current = gfIDs_angio.get(i);

			int indx_mon = cmmFunc.searchListString_index(gfID_current,
					gfIDs_monocots);
			int indx_eud = cmmFunc.searchListString_index(gfID_current,
					gfIDs_eud);
			
//			List<List<String>> found = cmpRank.searchMap(mapGFAtha, gfID_current);
			
//			for(int k=0; k<found.size();k++){
//	
//				if(k==0){
//				
//				List<String> first_list = found.get(0);
//					
//				System.out.print(first_list.get(0) + "\t" + first_list.get(1) + "\t"
//				+ first_list.get(2).split(";")[0]+"\t");
			
			System.out.print(gfID_current+"\t");

			    System.out.print(ranks_angio.get(i) + "\t"
					+ ranks_eud.get(indx_eud) + "\t"
					+ ranks_monocots.get(indx_mon) + "\t"
					+ optR_angio.get(i) + "\t" + optR_eud.get(indx_eud)
					+ "\t" + optR_monocots.get(indx_mon) + "\t"
					+ optLambdas_angio.get(i) + "\t"
					+ optLambdas_eud.get(indx_eud)+"\t"
					+optLambdas_monocots.get(indx_mon) + "\t"
					+ opt_likelihoods_angio.get(i) + "\t"
					+opt_likelihoods_eud.get(indx_eud)+"\t"
					+ opt_likelihoods_monocots.get(indx_mon) + "\t"
					+ pvalues_angio.get(i) + "\t"
					+pvalues_eud.get(indx_eud)+"\t"
					+ pvalues_monocots.get(indx_mon));

			System.out.println();
				}
				
				
//				else{ //to print rest of Atha genes in next lines
//					List<String> list = found.get(k);
//					
//					System.out.print(list.get(0) + "\t" + list.get(1) + "\t"
//					+ list.get(2).split(";")[0]);
//					
//					System.out.println();		
//				}
//				}

//			}
		}
		
		

		
		
		/** for final combined output file*/
		



}
