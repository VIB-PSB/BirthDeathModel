
package be.ugent.psb.setas.bdmodel.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.TreeMap;

import be.ugent.psb.setas.bdmodel.parsers.NewickParser;
import be.ugent.psb.setas.bdmodel.parsers.ReadComOutput_optValues;
import be.ugent.psb.setas.bdmodel.parsers.ReadGFcountsFile;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;
import be.ugent.psb.setas.bdmodel.parsers.WGMparser;
import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class Pvalues {

	/**
	 * p-values are the number of generated Observations that have less
	 * likelihood than the reference likelihood + half of the ones with equal
	 * likelihood divided by the whole number of generated observations
	 */

	private Node root;
	private ArrayList<Node> aryln;
	private int numOfObservations;
	private GenerateObservations go;
	private LikeLihood lc;
	private double lambda;

	static final double EPSILON = 1e-60;
//	static final double EPSILON = 1e-5;

	public Pvalues(Node root, ArrayList<Node> aryln , double lambda, int numberOfObservations){
		
		System.out.println("pvalues is not using any Cache");
		this.root = root;
		this.aryln = aryln;
		this.numOfObservations = numberOfObservations;
		this.go = new GenerateObservations(0,100,false);
		this.lc = new LikeLihood(lambda,aryln.get(0).getmaxNodeSize()+1);	
	}
	
//	public Pvalues(Node root, ArrayList<Node> aryln , double lambda, int numberOfObservations, ProbCalculator probCache, int lMCMC ){
//		this.root = root;
//		this.aryln = aryln;
//		this.lambda= lambda;
//		this.numOfObservations = numberOfObservations;
//		this.go = new GenerateObservations(5,100,false,probCache,lMCMC);
//		this.lc = new LikeLihood(lambda,aryln.get(0).getmaxNodeSize()+1, probCache);	
//	}
//	

	public static boolean isEqual(double a, double b, double PRECISION){
		return (Math.abs(a-b) < PRECISION);
	}
	
	
	public double calculateConditionalPvalues(int s, int rootSize, double refLogLk) {

		/* The values s, of the main root size, i.e. angiosperm size is neede to refer to the correct slot in the likelihood array, 
		 * because the method likelihood calculates likelihoods for all root sizes anyway
		 * so we need this to know for which root size are we actually calculating the p-values
		 * */
		
		    double pValue =0;
		
			for (int randObs = 1; randObs <= numOfObservations; randObs++) {

				/* because the root node that we use to start generating observations recursively on the tree is "eudicots"
				 * */
//				rootSize = eudicotSize;
				
				int[] genObsArray = go.generateObservation(root, rootSize, this.lambda);
				
				root.setLeafValues(genObsArray);

				double[] lkOfThisObs = lc.calcInternalLk(aryln);	
				
				double [] logLkOfThisObs = MathOperations.giveLogArray(lkOfThisObs);
				
				if(isEqual(logLkOfThisObs[s],refLogLk,EPSILON)){
					
					pValue += 0.5;			
				}
				
				/* Because the obsrevations are genrated according to the random BD be.ugent.psb.setas.bdmodel.model, so the likelihoods should be more or equal, if not, it means that the current observatio
				 * is following a special pattern rather the one described by the be.ugent.psb.setas.bdmodel.model, so in this case we penalize the p-value and increase it*/
				else if (logLkOfThisObs[s] < refLogLk) { // if it was not equal according to the precision
					pValue += 1;		
				}
				
			}
//			System.out.println(pValue);
			pValue = pValue*1.000000/numOfObservations*1.000000;
		
		return pValue;
	}

//	public double[] calculatePvalues(double lambda, double [] refLk) {
//
//		LikeLihood lc = new LikeLihood();
//
//		double[] pValue = new double[refLk.length];
//		pValue[0] = 0;
//
//		for (int s = 1; s < pValue.length; s++) {
//
//			for (int r = 1; r <= numberOfObservations; r++) {
//
//				int[] genObsArray = go.generateObservation(root, s, lambda);
//
//				root.setLeafValues(genObsArray);
//
//				double[] lkOfThisObs = lc.calcInternalLk(aryln, lambda);
//
//				if (lkOfThisObs[s] < refLk[s]) {
//					pValue[s] += 1;
//				}
//
//				if (lkOfThisObs[s] == refLk[s]) {
//					pValue[s] += 0.5;
//				}
//
//			}
//
//			pValue[s] = (Double) (pValue[s] / numberOfObservations);
//		}
//
//		return pValue;
//	}

//	 public static void main(String[]args){
//		 
//		 // Careful: We have changed genSize class to compare alpha to 1: wrong move! 
//	
////	long startTime = System.currentTimeMillis();
//		CommonFunctions cmf = new CommonFunctions();
//			
//		int maxNodeSize = 100;
//		double partitionSize = 0.1;
//			
//		//args: 0= tree, 1= wgd file,
//		Node root = SpeciesTreeParser.buildInsertWGDsandPartitionTree(args[0], args[1],
//		partitionSize, maxNodeSize);
//		ArrayList<Node> leaves = root.getLeaves();
//		
//		Queue<Node> queue = root.postOrder();
//		ArrayList<Node> aryln = new ArrayList<Node>(queue);
//
////		ReadComOutput_optValues rdopt = new ReadComOutput_optValues();
////		List<List<Double>> table = new ArrayList<List<Double>>();	
////		table = rdopt.readInputFile(args[2]);
////		ArrayList<String> gfIDs = rdopt.getGfIDs();
//
//		int numOfObs =1000;
//		String combinedOutput = args[2];
//		List<List<String>> map = cmf.readMapFile(combinedOutput);	
//		ProbCalculator probCalc = new ProbCalculator();
//		
//		int lMCMC =1000;
//		
//		int startGFindex = Integer.parseInt(args[3]);
//		int plusGF = Integer.parseInt(args[4]);
//		
//		for (int i = startGFindex; i <= startGFindex+plusGF; i++) {
//
//			List<String> gfrecord = map.get(i);
//			String gfID = gfrecord.get(0);
//
//			int rootSizeStar = Integer.parseInt(gfrecord.get(1));
//			double lambdaStar = Double.parseDouble(gfrecord.get(2));
//			double loglkStar = Double.parseDouble(gfrecord.get(3));
//			
//			 Pvalues pv= new Pvalues(root,aryln,lambdaStar,numOfObs, probCalc, lMCMC);
//			 double pValues = pv.calculateConditionalPvalues(rootSizeStar, rootSizeStar, loglkStar);
//			 System.out.print(gfID+"\t"+ lambdaStar+"\t"+pValues+"\n");
//			
//		}
//			

//		double wgtBlen = 0.260;
//		double oldrootBlen = 0.275776;
		
//		for(int i=0; i<table.size();i++){	
//		int i = Integer.parseInt(args[3]);
//
//		List<Double> record = table.get(i);	
//		double [] gf_optValues = new double [record.size()];
//		for (int j=0; j<gf_optValues.length;j++){
//			gf_optValues[j] = record.get(j);
//		}	
//		int r = (int)(Math.floor(gf_optValues[0]));
//		double lambda = gf_optValues[1];
//		double refLogLk = gf_optValues[2];
		
//		int eudicotSize = IncludeAngioSperms.generateSizeAtEudicots(root.getmaxNodeSize(), r, lambda,
//		wgtBlen, oldrootBlen);
		
//	    Pvalues pv= new Pvalues(root,aryln,lambda,numOfObs);
//        double pValues = pv.calculateConditionalPvalues(r, eudicotSize,
//		lambda, logLk);
//	    double pValues = pv.calculateConditionalPvalues(r, r, refLogLk);        
//		System.out.print(gfIDs.get(i)+"\t"+ pValues);
//		System.out.println("\n");
//		}

//		long endTime = System.currentTimeMillis();
//		System.out.print(endTime-startTime);
	
//	    double[] d= pv.calculateConditionalPvalues(s, eudicotSize, lambda, refLogLk);
//	 }

}
