package be.ugent.psb.setas.independent_parsers.KnKs;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class MakeNewKnKsTopBottomFileForMatlabPlots {

	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();

		// String matrixFile =
		// "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/KnKs/referenceLine/residuals/matrix_Zmay.txt";
		// String matrixFile =
		// "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/Exp_Func_newTopBottom1000/genefamily.angiosperm.core.Atha.codeml.Coex.TvsB";
		// String
		// matrixFile="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/Exp_Func_newTopBottom1000/TopBottom1000_ExpFunctPlots";
		String matrixFile = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/MusaFirst2WGDs/basedOnLambdaRanking/genefamily.angiosperm.core.txt.codeml";
		List<List<String>> map = cmf.readMapFile(matrixFile);

		// List<String> gfIDsInTheMapFile = cmf.readColX_String(matrixFile, 1);

		// String topGFidsFile =
		// "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/KnKs/repeat-KnKs_combRankLamBlock/fits_for_Lorenzo_20160509/input/TopBottomGFids/Top1000_combRankLamBlock_GFid.txt";
		// String bottomGFidsFile =
		// "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/KnKs/repeat-KnKs_combRankLamBlock/fits_for_Lorenzo_20160509/input/TopBottomGFids/Bottom1000_combRankLamBlock_GFid.txt";
		String topGFidsFile = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/MusaFirst2WGDs/basedOnCombinedRanking/TopBottom/Top1000GFids_combRank";
		String bottomGFidsFile = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/MusaFirst2WGDs/basedOnCombinedRanking/TopBottom/Bottom1000GFids_combRank";
		List<String> topGFids = cmf.read1ColFile_String(topGFidsFile);
		List<String> bottomGFids = cmf.read1ColFile_String(bottomGFidsFile);

		for (int i = 0; i < map.size(); i++) {// i=0; header

			List<String> mapRec = map.get(i);

			String mapGFid = mapRec.get(0);
			String speAbbreviation = mapRec.get(2).split("\\|")[0];

			if (topGFids.contains(mapGFid)) {

				for (int j = 0; j < mapRec.size(); j++) {
					System.out.print(mapRec.get(j) + "\t");
				}

				// System.out.print("\""+"Zmay"+"\""+"\t"+"\""+"top"+"\""+"\n");
				System.out.print(speAbbreviation + "\t" + "top" + "\n");

			}
			if (bottomGFids.contains(mapGFid)) {

				for (int j = 0; j < mapRec.size(); j++) {
					System.out.print(mapRec.get(j) + "\t");
				}

				// System.out.print("\""+"Zmay"+"\""+"\t"+"\""+"bottom"+"\""+"\n");
				System.out.print(speAbbreviation + "\t" + "bottom" + "\n");

			}
			
			if((!topGFids.contains(mapGFid)) && (!bottomGFids.contains(mapGFid))){
				
				
				for (int j = 0; j < mapRec.size(); j++) {
					System.out.print(mapRec.get(j) + "\t");
				}
				
				System.out.print(speAbbreviation + "\t" + "middle" + "\n");
				
			}
			//
			// to test:

			// List<String> keepTrack = new ArrayList<String>();

			// for(String top : topGFids){
			// for(String bottom: bottomGFids){
			//
			// if((!(gfIDsInTheMapFile.contains(bottom))) &&
			// (!(keepTrack.contains(bottom)))){
			//
			// keepTrack.add(bottom);
			// System.out.println("top gf not included: "+bottom);
			// }
			// if( (!topGFids.contains(mapGFid)) &&
			// (!bottomGFids.contains(mapGFid))){
			//
			// System.out.println(mapGFid);
			// }
			// }

		}

	}

}
