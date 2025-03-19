package be.ugent.psb.setas.independent_parsers.RVSEng;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class MakeASummaryFileRvsEngAllMGFs {

	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();
		String myRegex = "ORTHO";
		DecimalFormat df = new DecimalFormat("0.000000");
//		
//		int score =0;
//		String path = "/home/setas/Desktop/setas/Project1/RvsEng/BottomGFs/Tier1OutBottom_RvsEng_Rm1WGD/comOut_Bottom_rmWGD17";
//		String path = "/home/setas/Desktop/setas/Project1/RvsEng/MiddleGfs/RvsEng-MiddleGFs-37Angio9178";
		
		String path="/home/setas/Desktop/comOut_RvsEng28Eud_Top970";
		
		List<List<String>> map = cmf.readMapFile(path);
		int numberOfWGDs= 12;
		double [] loglkNew = new double[numberOfWGDs]; //for the 18 WGDs on the 37spe tree
		
		
		for (int j=1; j<=970; j++){
			System.out.println(j);
		}
		
		
//		for (int i = 0; i < map.size(); i += 3) {
		for (int i = 0; i < map.size(); i ++) {
				
			int score =0;
			double sumDiffLoglk=0;
			List<String> current = map.get(i);
			
			

//			if (current.get(0).split(myRegex).length > 1) {
//
//				String gfID = current.get(0);			
//				double rStar = Double.parseDouble(current.get(1));
//				double lambdaStar = Double.parseDouble(current.get(2));
//				double loglkStar = Double.parseDouble(current.get(3));
////				double pvStar = Double.parseDouble(current.get(4));
//								
//				System.out.print(gfID + "\t" + rStar + "\t" +df.format(lambdaStar) +"\t"
//						 + df.format(loglkStar) +"\t");
//						
//				for(int j=5; j< numberOfWGDs+5; j++){
//					
//					loglkNew [j-5] = Double.parseDouble(current.get(j));
//					
//					System.out.print(df.format(loglkNew[j-5])+"\t");
//					
//					if (loglkNew[j-5] < loglkStar){
//						
//						score+=1;
//					}
//					
//					sumDiffLoglk+= (loglkNew[j-5] - loglkStar);
//				}
//				
//				System.out.print(score+"\t");
////				double avgDiffLoglk = Double.parseDouble(current.get(25));
//				
//				double avgDiffLoglk= sumDiffLoglk / numberOfWGDs;
//				System.out.print(avgDiffLoglk+"\n");
//				
//				
////				List<String> fixedRlam = map.get(i + 1);
////				double logLkNew = Double.parseDouble(fixedRlam.get(3));
////				String result = fixedRlam.get(4);
//				
////				if (result.equals("lower")) {
////					score +=1;
////				}			
////				System.out.print(current.get(0) + "\t" + rStar + "\t" +df.format(lambdaStar) +"\t"
////						 + df.format(loglkStar) + "\t"+ df.format(logLkNew) + "\t"+ result +"\n");		
//			
//			}
			
	
			 
		}
	}

}
