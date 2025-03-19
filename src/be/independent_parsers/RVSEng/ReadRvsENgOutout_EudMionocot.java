package be.ugent.psb.setas.independent_parsers.RVSEng;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class ReadRvsENgOutout_EudMionocot {

	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();
		String path = "/home/setas/Desktop/setas/Project1/RvsEng/TopGFs/Monocots/Tier1Output/ComOut_EudMon_rmWGD5";

		List<List<String>> map = cmf.readMapFile(path);
		String myRegex = "ORTHO";

		List<List<String>> higherLogLowerLambda = new ArrayList<List<String>>();

		DecimalFormat df = new DecimalFormat("0.0000");

		/** To compare likelihood changes for a fixed rStar and lamStar */
		for (int i = 0; i < map.size(); i += 5) {

			List<String> current = map.get(i);

			if (current.get(0).split(myRegex).length > 1) {

				String gfID = current.get(0);
				double rStar = Double.parseDouble(current.get(1));
				double lambdaStar = Double.parseDouble(current.get(2));
				double loglkStar = Double.parseDouble(current.get(3));
				// double pvStar = Double.parseDouble(current.get(4));

				List<String> fixedRlam = map.get(i + 1);
				double newRefLogLk = Double.parseDouble(fixedRlam.get(3));

				List<String> fixedR_CPMlam = map.get(i + 2);
				double newLambdaMonocots = Double.parseDouble(fixedR_CPMlam.get(2));
				

				List<String> fixedRlam_Rm1WGD = map.get(i + 3);
				double loglk_Rm1WGD = Double.parseDouble(fixedRlam_Rm1WGD
						.get(3));
				String result = fixedRlam_Rm1WGD.get(4);

//				if (result.equals("lower")) {
//
//					System.out.print(gfID + "\t" + rStar + "\t"
//							+ df.format(lambdaStar) + "\t"
//							+ df.format(loglkStar) + "\t"
//							+ df.format(newRefLogLk) + "\t"
//							+ df.format(loglk_Rm1WGD) + "\t" + result);
//					System.out.println();
//				}

				
				List<String> CPM = map.get(i + 4);

				double RStar_CPM = Double.parseDouble(CPM.get(0));
				double lambdaStar_CPM = Double.parseDouble(CPM.get(1));
				double loglkStar_CPM = Double.parseDouble(CPM.get(2));
				
//				if(lambdaStar_CPM >= newLambdaMonocots){
//					
//					System.out.println(gfID+"\t"+rStar+"\t"+lambdaStar+"\t"+newLambdaMonocots+"\t"+RStar_CPM+"\t"+lambdaStar_CPM);
//				}
				
	           if(loglkStar_CPM <= newRefLogLk){
					
					System.out.println(gfID+"\t"+rStar+"\t"+loglkStar+"\t"+newRefLogLk+"\t"+RStar_CPM+"\t"+loglkStar_CPM);
				}

				// To check both logLk and lambda*New

				// List<String> newTreeRmWGD = map.get(i + 2);
				// double rNew = Double.parseDouble(newTreeRmWGD.get(0));
				// double lambdaNew = Double.parseDouble(newTreeRmWGD.get(1));
				// double logLkNew2 = Double.parseDouble(newTreeRmWGD.get(2));
				//
				// String compare ="lower";
				//
				// if(logLkNew2 > loglkStar){
				//
				// compare = "higher";
				// }
				// if (result.equals("lower") && lambdaNew > lambdaStar) {
				//
				// System.out.print(current.get(0) + "\t" + rStar + "\t" +
				// df.format(lambdaStar) +
				// "\t"
				// + df.format(loglkStar) + "\t"+ df.format(logLkNew) + "\t"
				// + result +
				// "\t"+rNew+"\t"+df.format(lambdaNew)+"\t"+df.format(logLkNew2)+"\t"+compare+"\n");
				//
				// }

			}
		}

		// /** To see the changes of lambda */
		// for (int i = 0; i < map.size(); i += 3) {
		//
		// List<String> current = map.get(i);
		//
		// if (current.get(0).split(myRegex).length > 1) {
		//
		// double rStar = Double.parseDouble(current.get(1));
		// double lambdaStar = Double.parseDouble(current.get(2));
		// double loglkStar = Double.parseDouble(current.get(3));
		// // double pvStar = Double.parseDouble(current.get(4));
		//
		// List<String> newTreeRmWGD = map.get(i + 2);
		//
		// double rNew = Double.parseDouble(newTreeRmWGD.get(0));
		// double lambdaNew = Double.parseDouble(newTreeRmWGD.get(1));
		// double logLkNew = Double.parseDouble(newTreeRmWGD.get(2));
		//
		// // if (lambdaNew > lambdaStar) {
		//
		// System.out.print(current.get(0) + "\t" + rStar + "\t"
		// + rNew + "\t" + lambdaStar + "\t" + lambdaNew
		// + "\t" + loglkStar + "\t" + logLkNew + "\n");
		// // }
		//
		// }
		// }

	}
}
