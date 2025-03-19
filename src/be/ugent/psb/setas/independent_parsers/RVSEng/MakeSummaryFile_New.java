package be.ugent.psb.setas.independent_parsers.RVSEng;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class MakeSummaryFile_New {

	// to compare the absolute values of the Loglks instead of the LRT test.
	public int countNumberOfNegativeLTRs(List<List<String>> map, int numOfGF, int numOfFieldsBeforeFirstNi, int numOfWGDs) {

		List<String> gf_record = map.get(numOfGF);
		String gfID = gf_record.get(0);
		
		System.out.print(gfID+"\t");
		int numOfWGDsWithLessLKAfterRm = 0;

		for (int j =0; j < numOfWGDs; j++) {
			
			double lrt = Double.parseDouble(gf_record.get(j+numOfFieldsBeforeFirstNi));

			if(lrt < 0){ /** if testing H0=full tree, the number of retrieved WGDs will be obtained when this is negative */
//			if (lrt > 0) { /** if H0 = rmWGD */
//			if (lrt - 2.00 < 0){ /** for AIC= Q <2  **/
			numOfWGDsWithLessLKAfterRm += 1;
			
			System.out.print("\t"+j);
			}
		}
		System.out.print("\n");
				
		return numOfWGDsWithLessLKAfterRm;
	}

	public int printSumOfLRTsBiggerThanThreshold(String pythonOutput_Ni_s, double threshold, int GFnumber, int numOfFieldsBeforeFirstNi) {

		CommonFunctions cmf = new CommonFunctions();
		List<List<String>> map = cmf.readMapFile(pythonOutput_Ni_s);

		int sum = 0;
		List<String> gf_record = map.get(GFnumber);

		int numberOfWGDs = 20;
		double[] n_i = new double[numberOfWGDs];

		int numOfWGDsInGFrecord = gf_record.size() - numOfFieldsBeforeFirstNi;

		for (int j = 0; j < numOfWGDsInGFrecord; j++) { /** to avoid errors for incomplete records **/

			n_i[j] = Double.parseDouble(gf_record.get(j + numOfFieldsBeforeFirstNi)) / 1000; /** 1=rootsize and 2=lambda **/

			if (n_i[j] > threshold) { /** n_i[j] > 0.95 , we can accept the null hypothesis of full-Tree **/

				sum += 1;
			}
		}

		return sum;

	}
	
	
//	public int countGFsThatRetrieveAWGD(int numOfColumn, String mapFile){
//		
//		
//		
//	}

	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();
		MakeSummaryFile_New mksum = new MakeSummaryFile_New();
		int numOfFieldsBeforeFirstNi = 3; //for full tree file: gfID, rootSize, lambda
//		int numOfFieldsBeforeFirstNi =2; //for rmWGD: gfID , rootsize
		int numOfWGDs =20;

		// String mapFile =
		// "/home/setas/Desktop/setas/RvsEng/Musa3TauMon2/H0isFullTree/python_output_ni";
//		String mapFile = "/home/setas/Desktop/setas/RvsEng/Musa3TauMon2/H0isRmWGD/python_ni_rmWGD";
//		List<List<String>> map = cmf.readMapFile(mapFile);
		
		String originalLRTfile ="/home/setas/Desktop/setas/RvsEng/Musa3TauMon2/H0isFullTree/OriginalLRTs";
//		String originalLRTfile ="/home/setas/Desktop/setas/RvsEng/Musa3TauMon2/H0isRmWGD/originalLRTs_merged_rmWGDs";
		List<List<String>> originalLRTsmap = cmf.readMapFile(originalLRTfile);

		//per family
//		for (int i = 0; i < 9178; i++) {
			//
//			List<String> gf_record = originalLRTsmap.get(i);
//			String gfID = gf_record.get(0);
			//
			// int numRetrievedWGDs =
			// mksum.printSumOfLRTsBiggerThanThreshold(mapFile, 0.95, i);
			// System.out.println(gfID + "\t" + numRetrievedWGDs);

//			int numOfWGDsWithLessLkAfterRm = mksum.countNumberOfNegativeLTRs(originalLRTsmap, i, numOfFieldsBeforeFirstNi, numOfWGDs);			
//			System.out.println(gfID+"\t"+numOfWGDsWithLessLkAfterRm);
			
			
			//
			// // if (sum_WGDs > 7) {
			// //
			// // System.out.print(gfID+"\t"+sum_WGDs);
			// //
			// // for (int j = 0; j < 20; j++) {
			// //
			// // if (n_i[j] > 0.95) {
			// // System.out.print("\t" + "wgd"+"\t"+j);
			// // }
			// // }
			// // System.out.print("\n");
			// // }
			//
			
//		}
			
			
		// per WGD
		
		for(int column =3; column <23 ; column++){
			
//		int column =3;
		int numberOfWGD = column -3; 
		List<Double> lksForWGDi = cmf.readColX_double(originalLRTfile, column);
		int sumNumberOfGFs =0;
		
//		int endLKinList = lksForWGDi.size();
		int endLKinList = 700;
			
		for(int k=0; k< endLKinList; k++){
						
			if(lksForWGDi.get(k) < 0){
				
				sumNumberOfGFs +=1;
			}
			
		}
		
		System.out.println(numberOfWGD +"\t"+ sumNumberOfGFs);
		}

		// for (int wgd=3; wgd<23; wgd++){
		//
		// ArrayList<Double> WGD = cmf.readColX_double(mapFile, wgd);
		// int sum_gfs = 0;
		//
		// for(double d : WGD){
		// if(d/1000 > 0.95){
		// sum_gfs +=1;
		// }
		// }
		//
		// System.out.println(sum_gfs);
		// }

	}

}
