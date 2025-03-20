package be.ugent.psb.setas.bdmodel.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import be.ugent.psb.setas.bdmodel.model.RVS_Engineering.CalculateCorrectLRToriginalFiles;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;
import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class GeneralErrorAssessment {

	Random rand = new Random();
	CalculateCorrectLRToriginalFiles calcLrt = new CalculateCorrectLRToriginalFiles();


	/** returns the index of next lambda which is > 0.1 + current lambda */
	public int returnLastIndex(int startIndex, double[] lambdas) {

		int endIndex = startIndex + 1;

		while ((endIndex < lambdas.length) && lambdas[endIndex] < lambdas[startIndex] + 0.1) {

			endIndex += 1;
		}
		return endIndex;
	}

	public int generateRandomCountToAdd_Normal(double mu, double std) {

		/** double r1 = mu+ rand.nextGaussian()*std; */
		double r1 = mu + rand.nextGaussian() * std;
		int r = (int) Math.floor(r1);

		if (r1 < 0) { // to have (-1,0] in the bin of zero, (-2,-1] in the bin
						// of -1, and so on..

			r += 1;
		}

		return r;
	}

	public int[] calcShuffledCounts(int[] originalObs, int minCount, int maxCount, double mu, double std) {

		int[] shuffeledObs = new int[originalObs.length];

		for (int i = 0; i < originalObs.length; i++) {

			int r = generateRandomCountToAdd_Normal(mu, std);
			shuffeledObs[i] = originalObs[i] + r;

			if (shuffeledObs[i] < minCount) {
				shuffeledObs[i] = minCount;
			}

			else if (shuffeledObs[i] > maxCount) {

				shuffeledObs[i] = maxCount;
			}
		}

		return shuffeledObs;
	}

	public double[] calc100LambdasPerInterval(int startIndex, double[] lambdasOrderedArray, List<String> GFidsOrdered,
			int[] rootSizes, String gfCountFile, int minCount, int maxCount, ProbCalculator probCache, Node root,
			int defMaxNodeSize, int numObsInEachBin, double mu, double std) {

		int endIndex = returnLastIndex(startIndex, lambdasOrderedArray);

		double[] newLambdas = new double[numObsInEachBin];

		for (int k = 0; k < numObsInEachBin; k++) {

			/* random = from + rndGenerator.nextInt(to - from + 1) */
			int randIndexGF = startIndex + rand.nextInt(endIndex - startIndex + 1);
			String gfID_chosen = GFidsOrdered.get(randIndexGF);
			int[] originalCounts = calcLrt.findObservationBasedOnGFid(gfCountFile, gfID_chosen, 37);
			int rootSize = rootSizes[randIndexGF];
			int[] shuffeledCOunts = calcShuffledCounts(originalCounts, minCount, maxCount, mu, std);
			
			int [] difference = new int[shuffeledCOunts.length];
			
			for(int i=0; i<shuffeledCOunts.length;i++){
			difference [i] = shuffeledCOunts[i]-originalCounts[i];
			}

			int [] hits = countHits(difference);
			int max  = returnAbsMax(difference);
			
			System.out.print(max);
			
			for(int h:hits){
				
				System.out.print("\t"+h);
			}
			
			System.out.print("\n");
			
//			System.out.println();
//			 for(int i=0; i< shuffeledCOunts.length;i++){
//			 System.out.print(shuffeledCOunts[i]-originalCounts[i]+"\t");
//			 }
//			 System.out.print("\n");
//			newLambdas[k] = calculateOptLambda(shuffeledCOunts, root, rootSize, defMaxNodeSize, probCache); // comment
																											// it
																											// out
																											// to
																											// only
																											// do
																											// analysis
																											// on
																											// the
																											// gene
																											// counts
		}

		return newLambdas;

	}

	public double calculateOptLambda(int[] shuffeledCounts, Node root, int rootSize, int defMaxNodeSize,
			ProbCalculator probCache) {

		double stepSize = 1e-4;// step calculating derivative
		double deltaLocalMoves = 1e-2;
		double tolD = 1e-3;
		double tolF = 1e-4;
		double minInterval = 1e-2;
		double maxInterval = 9.9999;
		double minAllowed = 1e-2;
		double maxAllowed = 9.9999;
		double precisionLambda = 1e-5; // one digit more than the number of
										// digits we require accurately
		int defaultmaxNodeSize = 100;
		root.setLeafValues(shuffeledCounts);
		ArrayList<Node> speciesTree = SpeciesTreeParser.setMaxNodeSize(root, defaultmaxNodeSize);

		CuttingPlaneMethod cpm = new CuttingPlaneMethod(speciesTree, rootSize, stepSize, deltaLocalMoves, tolD, tolF,
				minInterval, maxInterval, minAllowed, maxAllowed, precisionLambda, probCache);

		cpm.findOptimalLambda();
		return cpm.getOptimalLambda();
	}
	
	public int[] countHits(int [] shcounts){
		
		int [] countHits= new int[17];
		int[] refNumbers = new int[17];
		
		for(int i=0; i<17; i++){
			
			refNumbers[i] = i-8;
			
		}
		
		for(int s: shcounts){
			
			for (int r : refNumbers){
					
				if (s==r){
					
					countHits[s+8]+=1;
				}
				
			}
		}
		
		return countHits;
		
	}
	
	
	public int returnAbsMax(int [] myArray){
		
		int max= Math.abs(myArray[0]);
	
		for(int ma : myArray){
			
			if(Math.abs(ma) > Math.abs(max)){
				
				max = Math.abs(ma);
			}
		}
		
		return max;
	}

	public static void main(String[] args) {

		GeneralErrorAssessment genErr = new GeneralErrorAssessment();
		CommonFunctions cmf = new CommonFunctions();
		ProbCalculator probCache = new ProbCalculator();

		int minCount = 0;
		int maxCount = 100;
		double partitionSize = 0.1;
		int defaultMaxNodeSize = 100;
		int numObsInEachBin = 100;
		
		double mu=0;
		double std= 0.50;

		 String gfCountFile = args[2];

		Node root = SpeciesTreeParser.buildInsertWGDsandPartitionTree(args[0], args[1], partitionSize, defaultMaxNodeSize);
		ArrayList<Node> leaves = root.getLeaves();

		String combinedOutputFile_ordered = args[3];
		List<String> GFidsOrdered = cmf.readColX_String(combinedOutputFile_ordered, 0);
		List<Integer> rootSizes = cmf.readColX_Int(combinedOutputFile_ordered, 1);
		List<Double> lambdasOrdered = cmf.readColX_double(combinedOutputFile_ordered, 2); /**combinedOutput text file must be ordered on the column of lambdas */


		double[] lambdasOrderedArray = new double[lambdasOrdered.size()];
		for (int i = 0; i < lambdasOrdered.size(); i++) {

			lambdasOrderedArray[i] = lambdasOrdered.get(i).doubleValue();
		}

		int[] rootSizesArray = new int[rootSizes.size()];
		for (int j = 0; j < rootSizes.size(); j++) {
			rootSizesArray[j] = rootSizes.get(j).intValue();
		}

		int startIndex = Integer.parseInt(args[4]);

		while (startIndex < 9177) {

			int endIndex = genErr.returnLastIndex(startIndex, lambdasOrderedArray);

//			System.out.println(startIndex + "\t" + endIndex);

			double[] newLambdas = genErr.calc100LambdasPerInterval(startIndex, lambdasOrderedArray, GFidsOrdered,
					rootSizesArray, gfCountFile, minCount, maxCount, probCache, root, defaultMaxNodeSize,
					numObsInEachBin,mu,std);

//			for (double lam : newLambdas) {
//				System.out.println(lam);
//			}

			startIndex = endIndex + 1;
		}

	}

}
