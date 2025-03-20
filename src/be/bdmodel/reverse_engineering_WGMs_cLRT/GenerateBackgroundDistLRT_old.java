package be.ugent.psb.setas.bdmodel.model.RVS_Engineering;

import java.util.ArrayList;
import be.ugent.psb.setas.bdmodel.model.LikeLihood;
import be.ugent.psb.setas.bdmodel.model.MathOperations;
import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.bdmodel.model.ProbCalculator;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;

public class GenerateBackgroundDistLRT_old {

	private double lambdaPartitionSize;
	private int rootSize;
	private ProbCalculator probCalc;

	public GenerateBackgroundDistLRT_old(){}
	
	public GenerateBackgroundDistLRT_old(ProbCalculator pCalc) {

		this.probCalc = pCalc;}
	
	public GenerateBackgroundDistLRT_old(double lambdaPartitionSize, int r, ProbCalculator pCalc) {

		this.lambdaPartitionSize = lambdaPartitionSize;
		this.rootSize = r;
		this.probCalc = pCalc;
	}

	public double[] calcLambdasForSim(double minLamb, double maxLam) {

		int numOfPartitions = ((int) ((maxLam - minLamb) / lambdaPartitionSize))+2;
		double[] lambdas = new double[numOfPartitions];
		lambdas[0] = minLamb + (lambdaPartitionSize / 2);

		for (int i = 1; i < lambdas.length; i++) {

			lambdas[i] = lambdas[0] + i * (lambdaPartitionSize);}
		return lambdas;
	}

	public double calcLoglk(Node root, double lambda, int[] observation) {

		root.setLeafValues(observation);
		ArrayList<Node> aryln = root.postOrder(root);
		LikeLihood lk = new LikeLihood(lambda, 100, probCalc);

		double[] lks = lk.calcInternalLk(aryln);
		double[] logLks = MathOperations.giveLogArray(lks);
		return (logLks[rootSize]);
	}

	public ArrayList<Node> rootRmWGDs(String tree, String WGDall, String WGD8, String WGD9, String WGD12, String WGD16,
			int numberWGDs, double partitionSizeBranch, int defaultmaxNodeSize) {

		ArrayList<Node> rootRmWGDs = new ArrayList<Node>();

		for (int ignoreLineNumWGD = 0; ignoreLineNumWGD < numberWGDs; ignoreLineNumWGD++) {

			String WGDfile = WGDall;
			if (ignoreLineNumWGD == 8) {WGDfile = WGD8;}
			if (ignoreLineNumWGD == 9) {WGDfile = WGD9;}
			if (ignoreLineNumWGD == 12) {WGDfile = WGD12;}
			if (ignoreLineNumWGD == 16) {WGDfile = WGD16;}

			Node rRmWGD = SpeciesTreeParser.buildAndPartitionTree_ReverseEng(tree, WGDfile, partitionSizeBranch,
					defaultmaxNodeSize, ignoreLineNumWGD);
			rootRmWGDs.add(ignoreLineNumWGD, rRmWGD);
		}
		return rootRmWGDs;
	}

//	public static void main(String[] args) {
//
//		// args: 0=tree , 1 = WGDs(normal), 2 = Rm8, 3= Rm9, 4=Rm12 , 5=Rm16, 6
//		// = LambdaIndex
//		int rootSize = 3;
//		int defaultmaxNodeSize = 100;
//		double partitionSize = 0.1;
//		int lengthMCMC = 1000;
//		int numOfObservations = 1000;
//		int numberWGDs = 20;
//
////		double minLam = 0.01;
////		double maxLam = 0.10;
//		double partitionSizeLambda = 0.1;
//
//		ProbCalculator probCalc = new ProbCalculator();
//		GenerateBackgroundDistLRT genBgLRT = new GenerateBackgroundDistLRT(partitionSizeLambda, rootSize, probCalc);
//
//		Node rootOriginal = SpeciesTreeParser.buildInsertWGDsandPartitionTree(args[0], args[1], partitionSize,
//				defaultmaxNodeSize);
//
//		ArrayList<Node> rootRmWGDs = genBgLRT.rootRmWGDs(args[0], args[1], args[2], args[3], args[4], args[5],
//				numberWGDs, partitionSize, defaultmaxNodeSize);
//
////		double[] lambdas = genBgLRT.calcLambdasForSim(minLam, maxLam);
////		int index = Integer.parseInt(args[6]); // from index of lambda-0 to
//												// index of lambda-max, eg. for
//												// R=1 index= 0-16 (PBSarrayID=
//												// 3-19 to reflect the real
//												// values of lambdas according
//												// to minimum lambda =0.3) can
//												// loop on this for array jobs 
//												//
////		double lambdaCurrent = lambdas[index];
//		double lambdaCurrent = 4.441;
//
//		double[] loglkRmWGDs = new double[numberWGDs];
//		double[] lrtRmWGDs = new double[numberWGDs];
//		GenerateObservations go = new GenerateObservations(0, 100, false, probCalc, lengthMCMC);
//
//		for (int j = 0; j < numOfObservations; j++) {
//			
//			/** Ho: full-tree, H1: remove 1 WGD*/
//
////			int[] observation_H0 = go.generateObservation(rootOriginal, rootSize, lambdaCurrent);
////
////			double loglkFullTree = genBgLRT.calcLoglk(rootOriginal, lambdaCurrent, observation_H0);
////
////			for (int i = 0; i < numberWGDs; i++) {
////
////				loglkRmWGDs[i] = genBgLRT.calcLoglk(rootRmWGDs.get(i), lambdaCurrent, observation_H0);
////				lrtRmWGDs[i] = 2 * (loglkRmWGDs[i] - loglkFullTree);
////
////				System.out.print(lrtRmWGDs[i] + "\t");
////			}
////			System.out.print("\n");
//					
//			/** Ho: remove 1 WGD, H1: full-tree*/
//			for(int i=0; i<numberWGDs;i++){
//				
//				int[] observation_H0 = go.generateObservation(rootRmWGDs.get(i), rootSize, lambdaCurrent);
//
//				loglkRmWGDs[i] = genBgLRT.calcLoglk(rootRmWGDs.get(i), lambdaCurrent, observation_H0);
//				
//				double loglkFullTree = genBgLRT.calcLoglk(rootOriginal, lambdaCurrent, observation_H0);
//				
//				lrtRmWGDs[i] = 2 * (loglkFullTree -loglkRmWGDs[i]);
//				
//				System.out.print(lrtRmWGDs[i]+"\t");
//				
//			}		
//			System.out.print("\n");
//
//		}
//
//	}

}
