package be.ugent.psb.setas.independent_parsers.TandemBlock;


import java.util.ArrayList;

import be.ugent.psb.setas.bdmodel.model.MathOperations;
import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class TopBottomMixDataForBoxPlot {

	public static void main(String[] args) {

		CommonFunctions cmmFunct = new CommonFunctions();

		String pathBottom = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/TBdupPercent/onlyNonZeroFamilies/TopBottom/PT_Bottom";
		String pathTop = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/TBdupPercent/onlyNonZeroFamilies/TopBottom/PT_Top";

		ArrayList<Double> isBlock_B = cmmFunct.readColX_double(pathBottom, 1);
		ArrayList<Double> isTandem_B = cmmFunct.readColX_double(pathBottom, 2);
		ArrayList<Double> isBT_B = cmmFunct.readColX_double(pathBottom, 3);

		int lenxBottom = isBlock_B.size();

		ArrayList<Double> isBlock_T = cmmFunct.readColX_double(pathTop, 1);
		ArrayList<Double> isTandem_T = cmmFunct.readColX_double(pathTop, 2);
		ArrayList<Double> isBT_T = cmmFunct.readColX_double(pathTop, 3);

		int lenxTop = isBlock_T.size();
		
		double meanBlockTop = MathOperations.calcAverage(isBlock_T);
		double meanTandemTop=MathOperations.calcAverage(isTandem_T);
		double meanBTtop = MathOperations.calcAverage(isBT_T);
		
		double stdBlockTop = MathOperations.calcStd(isBlock_T);
		double stdTandemTop = MathOperations.calcStd(isTandem_T);
		double stdBTtop = MathOperations.calcStd(isBT_T);

		
		double meanBlockBottom = MathOperations.calcAverage(isBlock_B);
		double meanTandemBottom=MathOperations.calcAverage(isTandem_B);
		double meanBTbottom = MathOperations.calcAverage(isBT_B);
		
		double stdBlockBottom = MathOperations.calcStd(isBlock_B);
		double stdTandemBottom = MathOperations.calcStd(isTandem_B);
		double stdBTbottom = MathOperations.calcStd(isBT_B);
		
		
		System.out.println(meanBlockTop+"\t"+meanBlockBottom);
		System.out.println(stdBlockTop+"\t"+stdBlockBottom);
		
		System.out.println(meanTandemTop+"\t"+meanTandemBottom);
		System.out.println(stdTandemTop+"\t"+stdTandemBottom);
		
		
		System.out.println(meanBTtop+"\t"+meanBTbottom);
		System.out.println(stdBlockTop+"\t"+stdBTbottom);
		
//		int lenx = Math.max(lenxTop, lenxBottom);

//		for (int i = 0; i < lenx; i++) {
//
//			if (i < lenxBottom && i < lenxTop) {
//
//				System.out.println(isBlock_T.get(i) + "\t" + isBlock_B.get(i)
//						+ "\t" + isTandem_T.get(i) + "\t" + isTandem_B.get(i)
//						+ "\t" + isBT_T.get(i) + "\t" + isBT_B.get(i));
//			}
//
//			if (i < lenxBottom && i >= lenxTop) {
//
//				System.out.println("" + "\t" + isBlock_B.get(i)
//						+ "\t" + "" + "\t" + isTandem_B.get(i)
//						+ "\t" + "" + "\t" + isBT_B.get(i));
//			}
//
//			if (i >= lenxBottom && i < lenxTop) {
//
//				System.out.println(isBlock_T.get(i) + "\t" + ""
//						+ "\t" + isTandem_T.get(i) + "\t" + ""
//						+ "\t" + isBT_T.get(i) + "\t" + "");
//			}
//
//		}

	}

}
