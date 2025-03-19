package be.ugent.psb.setas.bdmodel.oldclasses;

import java.util.ArrayList;
import java.util.Random;

import be.ugent.psb.setas.bdmodel.model.LikeLihood;
import be.ugent.psb.setas.bdmodel.model.Node;

/**
 * Simulated Annealing to find the optimal lambda for a given root size and
 * observation array.
 */
public class SimulatedAnnealing {

	public Node root;
	public double tempZero = 1000;
	public int numOfSteps = 1030;
	public double beta = 0.99;
	double stepSize = 0.01;
	Random random = new Random();

	/*
	 * assume we fix max (s+c) = 200 i.e. s<101 and c<101. (maxNodeSize=100 but lk
	 * array size will be doubled) for alphas <= 0.520336 we have no probability out
	 * of range. which means: lambda*t <= 1.084792688 .
	 */

	// double randomValue = rangeMin + (rangeMax - rangeMin) * r.nextDouble();
	// testLambda = (1.084792688 / maxBranchLength) * random.nextDouble();
	// testLambda = (1.0 / maxBranchLength) * random.nextDouble();

	public double simulatedAnnealingForlambda(int rootSize, double lambdaZero) {

		ArrayList<Node> arl = root.postOrder(root);

		// double maxBranchLength = root.maxBlen();
		double[] lambdaSeq = new double[numOfSteps];

		double[] temperature = new double[numOfSteps];
		lambdaSeq[0] = lambdaZero;
		temperature[0] = tempZero;

		double convergingValue = lambdaZero;

		LikeLihood lc = new LikeLihood(lambdaZero, root.getmaxNodeSize() + 1);
		double[] oldLk = lc.calcInternalLk(arl);
		double oldLik = oldLk[rootSize];

//	 int counter = 0;

		for (int t = 1; t < numOfSteps; t++) {
//	 if (counter >= 1000) {
//	 break;
//	 } 
//	 System.out.println("t: "+t);
//     System.out.println("old likelihood: "+ oldLik);

			temperature[t] = tempZero * Math.pow(beta, t);
			double u = random.nextDouble();

//	 double min = lambdaSeq[t-1]+Math.log10((temperature[t]-1)/10000);
//	 
//	 if( min <0 ){
//		 min = 0;
//	 }
//	 double max = lambdaSeq[t-1]-Math.log10((temperature[t]-1)/10000); 
//	 if (max > 1){		 
//		 max = 1;
//	 }
//	 double testLambda = min + (max-min)* random.nextDouble(); 

			double testLambda = lambdaZero + (t * stepSize);
//	 System.out.println("be.ugent.psb.setas.bdmodel.test lambda: "+ testLambda);

			LikeLihood lc_new = new LikeLihood(testLambda, root.getmaxNodeSize() + 1);
			double[] newlk = lc_new.calcInternalLk(arl);
			double newLik = newlk[rootSize];

//	 System.out.println("new likelihood: "+ newLik);

			double exponential = Math.exp((newLik - oldLik) / temperature[t]);

			double alpha = Math.min(1, exponential);

			if (alpha >= u) {
//	  System.out.println("accepted");  
				lambdaSeq[t] = testLambda;

//	 if (t > 5000 && (lambdaSequnece[t] - lambdaSequnece[t - 1]) <= 0.0001) {
//	 counter += 1;}

				/* update the new initial Likelihood equal to the new likelihood */
				oldLik = newLik;
			}

			if (alpha < u) {

//	 System.out.println("rejected");
				lambdaSeq[t] = lambdaSeq[t - 1];
				/* keep the new initial Likelihood equal to the old likelihood */}
			convergingValue = lambdaSeq[t];
		}
		return convergingValue;
	}

	/*
	 * Becasue it may be fatser to store and use all optimal lambdas for different
	 * root sizes and try to compae all likelihoods directly --> brute force
	 * approach , to make sure if previous results still stay valid, but for now
	 * there's a bug in usage in the main methid
	 */

//	public double[] simulatedAnnealingForlambdaTwo(double lambdaZero) {
//
//		// double maxBranchLength = root.maxBlen();
//		int nodeSize = root.lk.length;
//		double[][] lambdaSequnece = new double[nodeSize][numOfSteps];
//
//		double[][] temperature = new double[nodeSize][numOfSteps];
//
//		lambdaSequnece[0][0] = lambdaZero;
//		temperature[0][0] = tempZero;
//
//		double testLambda = 0;
//		double[] convergingValue = new double[nodeSize];
//
//		for (int w = 1; w < nodeSize; w++) {
//
//			double[] oldLik = lc.calcInternalLk(arl, lambdaZero);
//
//			int counter = 0;
//
//			for (int t = 1; t < numOfSteps; t++) {
//
//				if (counter >= 1000) {
//					break;
//				}
//
//				temperature[w][t] = tempZero * Math.pow(beta, t);
//
//				double u = random.nextDouble();
//				testLambda = random.nextDouble();
//
//				double[] newLik = lc.calcInternalLk(arl, testLambda);
//
//				double exponential = Math.exp((newLik[w] - oldLik[w])
//						/ temperature[w][t]);
//
//				double alpha = Math.min(1, exponential);
//
//				if (alpha >= u) {
//					/*
//					 * update the new initial Likelihood equal to the new
//					 * likelihood
//					 */
//					lambdaSequnece[w][t] = testLambda;
//
//					if (t > 5000
//							&& (lambdaSequnece[w][t] - lambdaSequnece[w][t - 1]) <= 0.001) {
//						counter += 1;
//					}
//					oldLik = newLik;
//				}
//
//				if (alpha < u) {
//					lambdaSequnece[w][t] = lambdaSequnece[w][t - 1];
//					/*
//					 * keep the new initial Likelihood equal to the old
//					 * likelihood
//					 */
//				}
//				convergingValue[w] = lambdaSequnece[w][t];
//			}
//
//		}
//		return convergingValue;
//	}

//	 public static void main(String[] args) {
//	
//		 
//		 
//		 NewickParser np = new NewickParser();
//			Node root = np
//					.buildTree("/home/setas/workspace/BirthDeathModel/src/Files/EudicotTree-Arab");		
//			 WGMparser wgm = new WGMparser();
//			 List<List<String>> wgdList =
//			 wgm.readInputFile("/home/setas/workspace/BirthDeathModel/src/Files/WGDeudicotTree-Arab");
//			 root.insertWGM(root, wgdList);					
//			ArrayList<Node> leaves = new ArrayList<Node>();
//			leaves = root.getLeaves();
//			root.setNumberOfLeaves(leaves.size());
//	
//	 int[] idealObs = new int[10];
//		 idealObs[0] = 10;
//		 idealObs[1] = 20;
//		 idealObs[2] = 10;
//		 idealObs[3] = 6;
//		 idealObs[4] = 8;
//		 idealObs[5] = 5;
//		 idealObs[6] = 10;
//		 idealObs[7] = 5;
//		 idealObs[8] = 5;
//		 idealObs[9] = 5;
//	
//	 root.setLeafValues(idealObs);
//	
//	 SimulatedAnnealing sa = new SimulatedAnnealing();
//	 sa.root = root;
//	 double lambdaZero = 0;
//	
//	 double optimalLambda =
//	 sa.simulatedAnnealingForlambda(5, lambdaZero);
//	 System.out.println("optimal lambda is: "+optimalLambda);
//	
//	 }
}
