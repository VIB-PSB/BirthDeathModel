package be.ugent.psb.setas.independent_parsers.RVSEng;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class ReadReverseEngOutput {

	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();
		String path = "/home/setas/Desktop/setas/Project1/RvsEng/TopGFs/AngioSperms/individualMGFs/Tier1Output_RvsEng_Rm1WGD/ComOut_WGD0";

		List<List<String>> map = cmf.readMapFile(path);
		String myRegex = "ORTHO";

		List<List<String>> higherLogLowerLambda = new ArrayList<List<String>>();

		DecimalFormat df = new DecimalFormat("0.000");

		/** To compare likelihood changes for a fixed rStar and lamStar */
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
		// List<String> fixedRlam = map.get(i + 1);
		// double logLkNew = Double.parseDouble(fixedRlam.get(3));
		// String result = fixedRlam.get(4);
		//
		// // if (result.equals("lower")) {
		// //
		// // System.out.print(current.get(0) + "\t" + rStar + "\t" +
		// // df.format(lambdaStar) +"\t"
		// // + df.format(loglkStar) + "\t"+ df.format(logLkNew) + "\t"+ result
		// +"\n");
		// //
		// // }
		//
		// // To check both logLk and lambda*New
		//
		// // List<String> newTreeRmWGD = map.get(i + 2);
		// // double rNew = Double.parseDouble(newTreeRmWGD.get(0));
		// // double lambdaNew = Double.parseDouble(newTreeRmWGD.get(1));
		// // double logLkNew2 = Double.parseDouble(newTreeRmWGD.get(2));
		// //
		// // String compare ="lower";
		// //
		// // if(logLkNew2 > loglkStar){
		// //
		// // compare = "higher";
		// // }
		// // if (result.equals("lower") && lambdaNew > lambdaStar) {
		// //
		// // System.out.print(current.get(0) + "\t" + rStar + "\t" +
		// // df.format(lambdaStar) +
		// // "\t"
		// // + df.format(loglkStar) + "\t"+ df.format(logLkNew) + "\t"
		// // + result +
		// "\t"+rNew+"\t"+df.format(lambdaNew)+"\t"+df.format(logLkNew2)+"\t"+compare+"\n");
		// //
		// // }
		//
		// }
		// }

		/** To see the changes of lambda */
		for (int i = 0; i < map.size(); i += 3) {

			List<String> current = map.get(i);

			if (current.get(0).split(myRegex).length > 1) {

				String gfID = current.get(0);
				double rStar = Double.parseDouble(current.get(1));
				double lambdaStar = Double.parseDouble(current.get(2));
				double loglkStar = Double.parseDouble(current.get(3));
				// double pvStar = Double.parseDouble(current.get(4));

//				List<String> newTreeRmWGD = map.get(i + 2);
//				double rNew = Double.parseDouble(newTreeRmWGD.get(0));
//				double lambdaNew = Double.parseDouble(newTreeRmWGD.get(1));
//				 double logLkNew = Double.parseDouble(newTreeRmWGD.get(2));

//				if (lambdaNew > lambdaStar) {
//
//					System.out.print(gfID + "\t" + current.get(0) + "\t"
//							+ rStar + "\t" + rNew + "\t" + lambdaStar + "\t"
//							+ lambdaNew + "\t" + loglkStar + "\t" + logLkNew
//							+ "\n");
//				}
				// if (logLkNew < loglkStar) {
				//
				// System.out.print(gfID+"\t"+current.get(0) + "\t" + rStar +
				// "\t"
				// + rNew + "\t" + lambdaStar + "\t" + lambdaNew
				// + "\t" + loglkStar + "\t" + logLkNew + "\n");
				// }

				List<String> fixedRlam = map.get(i + 1);
				double logLkNew = Double.parseDouble(fixedRlam.get(3));
				String result = fixedRlam.get(4);
				 if (result.equals("lower")) {
				 System.out.print(current.get(0) + "\t" + rStar + "\t"
				 +df.format(lambdaStar) +"\t"
				 + df.format(loglkStar) + "\t"+ df.format(logLkNew) + "\t"+
				 result +"\n");
				
				 }

				// if (lambdaNew > lambdaStar && result.equals("lower")) {
				//
				// System.out.print(gfID+"\t"+current.get(0) + "\t" + rStar +
				// "\t"
				// + rNew + "\t" + lambdaStar + "\t" + lambdaNew
				// + "\t" + loglkStar + "\t" + logLkNew + "\n");
				// }

			}
		}

	}
}
