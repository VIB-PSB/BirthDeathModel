package be.ugent.psb.setas.bdmodel.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import be.ugent.psb.setas.bdmodel.model.RVS_Engineering.CalculateCorrectLRToriginalFiles;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;
import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class ErrorAssessment {
	
	public int generateNonZeroValueChange(int minVal, int maxVal){
		
		Random r = new Random();
		int valueChange = minVal+ r.nextInt(maxVal-minVal+1);
		
		if(valueChange==0){
			valueChange = minVal+ r.nextInt(maxVal-minVal+1);
		}
		
		return valueChange;
		
	}
	
	
	public int generateNonRepetetiveIndex (int numOfSpe, List<Integer> keepIndex){
		
		Random random = new Random();
		int randomIndex = random.nextInt(numOfSpe);
		
		 if(keepIndex.contains(randomIndex)){		 
			 randomIndex = random.nextInt(numOfSpe);
		 }	 
		return randomIndex;
	}

	public int[] generateObsWithError__valueFix_signMonotone(int[] originalgfCounts, double errorPercent, int valueChange,
			boolean isIncrement) {

		int numOfSpe = originalgfCounts.length;
		int [] shuffledObs = new int[originalgfCounts.length];
		System.arraycopy( originalgfCounts, 0, shuffledObs, 0, originalgfCounts.length );

		int numberOfcountsToChange = (int) (Math.floor(errorPercent * numOfSpe));

		// generate random indexes to change which corresponds to the random
		// species to change their counts
		/* random = from + rndGenerator.nextInt(to - from + 1) */
		List<Integer> keepRandIndex = new ArrayList<Integer>();

		for (int j = 0; j < numberOfcountsToChange; j++) {
		
			int randomIndex = generateNonRepetetiveIndex(numOfSpe,keepRandIndex);	

			if (isIncrement) {
				shuffledObs[randomIndex] = originalgfCounts[randomIndex] + valueChange;
			} else { // its decrement
				shuffledObs[randomIndex] = originalgfCounts[randomIndex] - valueChange;
			}
		}

		return shuffledObs;

	}

	public int[] generateObsWithError_valueFix_signRandom(int[] originalgfCounts, double errorPercent, int valueChange) {

		int numOfSpe = originalgfCounts.length;
		int [] shuffledObs = new int[originalgfCounts.length];
		System.arraycopy( originalgfCounts, 0, shuffledObs, 0, originalgfCounts.length );

		int numberOfcountsToChange = (int) (Math.floor(errorPercent * numOfSpe));

		Random random = new Random();
		List<Integer> keepRandIndex = new ArrayList<Integer>();

		for (int j = 0; j < numberOfcountsToChange; j++) {

			int randomIndex = generateNonRepetetiveIndex(numOfSpe,keepRandIndex);	

			int randomDecision = random.nextInt(10);

			if (randomDecision > 5) {
				shuffledObs[randomIndex] = originalgfCounts[randomIndex] + valueChange;
			}

			else {
				shuffledObs[randomIndex] = originalgfCounts[randomIndex] - valueChange;
			}

		}

		return shuffledObs;

	}
		
	public int[] generateObsWithError_valueRand_signRandom(int[] originalgfCounts, double errorPercent, int minVal, int maxVal, int minValid, int maxValid) {

		int numOfSpe = originalgfCounts.length;
		int [] shuffledObs = new int[numOfSpe];
		System.arraycopy( originalgfCounts, 0, shuffledObs, 0, numOfSpe );

		int numberOfcountsToChange = (int) (Math.floor(errorPercent * numOfSpe));

		Random random = new Random();
		List<Integer> keepRandIndex = new ArrayList<Integer>();
		
		for (int j = 0; j < numberOfcountsToChange; j++) {

			int randomIndex = generateNonRepetetiveIndex(numOfSpe,keepRandIndex);	
            int valueChange = generateNonZeroValueChange(minVal, maxVal);
			
			int randomDecisionSign = random.nextInt(10);

			if (randomDecisionSign > 5) {
				shuffledObs[randomIndex] = originalgfCounts[randomIndex] + valueChange;
			}

			else {
				shuffledObs[randomIndex] = originalgfCounts[randomIndex] - valueChange;
			}

		}

		return bringInValidRange(shuffledObs,minValid, maxValid);

	}
	
	public int[] bringInValidRange (int [] obs, int minValid, int maxValid){
		
		for(int i=0; i<obs.length;i++){
			
			if(obs[i] > maxValid){
				obs[i] = maxValid;
			}
			
			else if(obs[i] < minValid){
				obs[i] = minValid;
			}
		}
		
		return obs;
		
	}

	public static void main(String[] args) {

		CalculateCorrectLRToriginalFiles calcLrt = new CalculateCorrectLRToriginalFiles();
		ErrorAssessment err = new ErrorAssessment();
		CommonFunctions cmf = new CommonFunctions();

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
		double partitionSize = 0.1;
		int defaultmaxNodeSize = 100;
		
		int minValid =0;
		int maxValid = 100;

		Node root = SpeciesTreeParser.buildInsertWGDsandPartitionTree(args[0], args[1], partitionSize, defaultmaxNodeSize);
		ArrayList<Node> leaves = root.getLeaves();

		String combOutput = args[3]; // args3= comoutputfile
		List<List<String>> mapGFoptimalValues = cmf.readMapFile(combOutput);
		
		int gf= Integer.parseInt(args[4]);

		List<String> gfRecord = mapGFoptimalValues.get(gf);
		String gfID = gfRecord.get(0);
		int rootSizeStar = Integer.parseInt(gfRecord.get(1));
		double lambdaStar = Double.parseDouble(gfRecord.get(2));
		double lkStar = Double.parseDouble(gfRecord.get(3));
		System.out.println(gfID + "\t" + rootSizeStar);
		System.out.println(lambdaStar + "\t" + lkStar);

		int[] originalCounts = calcLrt.findObservationBasedOnGFid(args[2], gfID, 37);


		double errorPercent = Double.parseDouble(args[5]);		
//		int valueChange = 1;
//		boolean isIncrement = true;
		int minInterval_countErr= Integer.parseInt(args[6]);
		int maxInterval_countErr = Integer.parseInt(args[7]);
		int numberOfRepeats=1000;
		
		ProbCalculator probCache = new ProbCalculator();
		for(int repeat =0; repeat <numberOfRepeats; repeat++){

	    int[] shuffeledCounts = err.generateObsWithError_valueRand_signRandom(originalCounts, errorPercent, minInterval_countErr,maxInterval_countErr, minValid, maxValid);
	    
		root.setLeafValues(shuffeledCounts);
		ArrayList<Node> speciesTree = SpeciesTreeParser.setMaxNodeSize(root, defaultmaxNodeSize);

		 CuttingPlaneMethod cpm = new CuttingPlaneMethod(speciesTree,
		 rootSizeStar, stepSize, deltaLocalMoves, tolD,
		 tolF, minInterval, maxInterval, minAllowed, maxAllowed,
		 precisionLambda, probCache);
		
		 cpm.findOptimalLambda();
		
		 System.out.println(cpm.getOptimalLambda()+"\t"+cpm.getFoptimalLambda());

		 }
	}

}
