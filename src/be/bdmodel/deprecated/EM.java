package be.ugent.psb.setas.bdmodel.oldclasses;

import java.util.ArrayList;
import java.util.Queue;
import java.util.Random;
import java.util.Stack;

import be.ugent.psb.setas.bdmodel.model.FindMaxArray;
import be.ugent.psb.setas.bdmodel.model.LikeLihood;
import be.ugent.psb.setas.bdmodel.model.MathOperations;
import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.bdmodel.parsers.BuildTree;
import be.ugent.psb.setas.bdmodel.parsers.NewickParser;
import be.ugent.psb.setas.bdmodel.parsers.WGMparser;

public class EM {

	public Node root;
	double epsilon;

	FindMaxArray fm = new FindMaxArray();


	FindMaxArray fma = new FindMaxArray();


	
	public double emForLambda(int numberOfSteps) {
		
		ArrayList<Node> aryln = root.postOrder(root);
		SimulatedAnnealing sa = new SimulatedAnnealing();
		sa.root = root;

		/* To save results of each of 10 step of EM */
		double[] bestLambda = new double[numberOfSteps+1];
		bestLambda[0]=0;
		double[] bestLogLk = new double[numberOfSteps+1];
		bestLogLk[0]=0;
		int[] bestRootSize = new int[numberOfSteps+1];
		bestRootSize[0]=0;
		
		Random r = new Random();
//		double maxBranchLength = root.maxBlen();
		
		Stack<Double> stkLambda = new Stack<Double>();
		Stack<Double> stkLk = new Stack<Double>();
		Stack<Integer> stkMaxIndex = new Stack<Integer>();
		
		/* Initialization can not be arbitrary, since if in one oteration, we do not proceed 
		 * to any new lambda, the top of stack, which is the initial value would be returned,
		 * in other words, these values are added to the stack and could be st.peek()*/
		
		for (int i = 1; i < numberOfSteps+1; i++) {
			
			System.out.println("EM iteration:" + i);
			
			/* to start EM from different starting points */
			double lambdaZero = r.nextDouble();
			System.out.println("lambdaZero "+ lambdaZero);
			
			double nextEstimatedLambda = lambdaZero;
			stkLambda.add(nextEstimatedLambda);
				
			int counter =0;
			
			LikeLihood lk = new LikeLihood(lambdaZero, aryln.get(0).getmaxNodeSize()+1);
			double[] save = lk.calcInternalLk(aryln);			
			double [] logSave = MathOperations.giveLogArray(save); 

			double maxLogLk = fma.findMaxValue(logSave);
			int maxIndex = fma.findIndexOfMaxValue(logSave);
			stkMaxIndex.add(maxIndex);
			
			double initialLogLk = (-2)*epsilon + maxLogLk ;
			
			System.out.println("initial LogLikelihood "+ initialLogLk);
			stkLk.add(initialLogLk);
	
			while (((maxLogLk - stkLk.peek()) > epsilon) && (counter<100)) {

				System.out.println("next Estimated Lambda "+ nextEstimatedLambda);
				stkLambda.add(nextEstimatedLambda);
				
				System.out.println("max LogLikelihood "+ maxLogLk);
				stkLk.add(maxLogLk);
				
				System.out.println("max Index "+ maxIndex);
				stkMaxIndex.add(maxIndex);
					
				nextEstimatedLambda = sa.simulatedAnnealingForlambda(maxIndex,
						lambdaZero);
					
				LikeLihood lk_new = new LikeLihood(nextEstimatedLambda, aryln.get(0).getmaxNodeSize()+1);
				save = lk_new.calcInternalLk(aryln);
				logSave = MathOperations.giveLogArray(save); 
				
				maxLogLk = fma.findMaxValue(logSave);
				
				maxIndex = fma.findIndexOfMaxValue(logSave);
						
				counter +=1;
			}

			bestRootSize[i] = stkMaxIndex.peek();
			System.out.println("best root size at this step: "+ bestRootSize[i]);
			
			bestLambda[i] = stkLambda.peek();
			System.out.println("best Lambda at this step: "+ bestLambda[i]);
			
			bestLogLk[i] = stkLk.peek();
			System.out.println("best log lk at this step: "+ bestLogLk[i]);
			
//			nextEstimatedLambda = (1.084792688 / maxBranchLength)
//					* r.nextDouble();	
		}
	
		System.out.println("--------------------------------------------------------------------");
			
		int maxStep = fm.findIndexOfMaxValue(bestLogLk);
		System.out.println("optimal lk is obtained in step: "+maxStep+ " of EM");
		
		int optimalSize= bestRootSize[maxStep];
		System.out.println("optimal Size: "+optimalSize);
		
		root.setoptimalSize(optimalSize);
		
		double optimalLambda= bestLambda[maxStep];
		System.out.println("optimal Lambda: "+optimalLambda);
		
		return optimalLambda;
	}
}
	// public static void main(String[] args) {
	//
	// // // BuildTree bt = new BuildTree();
	// // // List<Node> list = bt.readInput(args[0]);
	// // // Node root = bt.buildTree(list, list.get(0));
	//
	// NewickParser np = new NewickParser();
	// Node root = np
	// .buildTree("/home/setas/workspace/BirthDeathModel/src/Files/plaza.txt.tree");
	//
	// // // WGD wgd = new WGD();
	// // // List<List<String>> wgdList = wgd.readInputFile(args[1]);
	// // // root.insertWgd(root, wgdList);
	//
	// ArrayList<Node> leaves = root.getLeaves(root);
	// // System.out.println("leafs: " + leaves.size());
	//
	// root.setNumberOfLeaves(leaves.size());
	// // System.out.println("nom of leaves: " + root.getNumberOfLeaves());
	//
	// int numberOfLeaves = root.getNumberOfLeaves();
	// int[] testObs = new int[numberOfLeaves];
	//
	// //  idealObs[0] = 4;
	// //  idealObs[1] = 4;
	// //  idealObs[2] = 2;
	// //  idealObs[3] = 4;
	//
	// for (int k = 0; k < numberOfLeaves; k++) {
	// testObs[k] = 1;
	// }
	// root.setLeafValues(root, testObs);
	//
	// //  ReadGFcountsFile rgf = new ReadGFcountsFile();
	// //  List<List<Integer>> counts = rgf.readInputFile(args[2]);
	// //
	// //  int numberOfGfs = counts.size();
	// //  to save the likelihood of each gene family counts:
	// //  double likelihoodOfGfs[] = new double[numberOfGfs];
	// //  int OptimalRootSizeOfaGF[] = new int [numberOfGfs];
	// //  double optimallambdaForGfandS[] = new double [numberOfGfs];
	// // // for (int i = 0; i < numberOfGfs; i++) {
	// //
	// //  List<Integer> li = counts.get(i);
	// //  root.setLeaveValues(List<Integer> valueList);
	// //  int[] originalObservation = new int[li.size()];
	// 
	// //  for (int m = 0; m < li.size(); m++) {
	// //  originalObservation[m] = li.get(m);
	// //  }
	//
	// EM em = new EM();
	//
	// int maxNodeSize=100;
	// root.setmaxNodeSize(maxNodeSize);
	// root.setMaxSizes(root, maxNodeSize);
	//
	// double lambdaZero = 0.001;
	// double epsilon = 0.0001;
	// double tempZero = 1000;
	// int numOfSteps = 10000;
	// double beta = 0.99;
	// //
	// int numberOfObservations = 100;
	// //
	// double bestLambda = em.emForLambda(lambdaZero);
	//
	// // int bestSize = em.emForRootSize(bestLambda, numberOfObservations);
	// // System.out.println("best root size: " + bestSize);
	// //
	// // LikeLihood lk = new LikeLihood();
	// // double[] optimizedLikelihoods = lk.calculateLikelihoods(root, 100,
	// // bestLambda, testObs);
	// // double optimizedLikelihood = optimizedLikelihoods[bestSize];
	// // System.out.println(optimizedLikelihood);
	// //
	// // // to sort the results with indexes:
	// // // //copy the array into a new array
	// // // double [] temp = likelihoodOfGfs;
	// // //
	// // // Arrays.sort(temp); //sort ascending
	// // //
	// // // //final array of indexes
	// // // int index_array[] = new int[numberOfGfs];
	// // //
	// // // //iteretate on x array
	// // // for(int i=0; i<numberOfGfs; i++)
	// // // //search the position of a value of the original x array into the
	// // // sorted y array, store the position in the index array
	// // // index_array[i] = Arrays.binarySearch(likelihoodOfGfs,temp[i]);
	//
	// }


