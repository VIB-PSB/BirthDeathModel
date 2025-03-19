
package be.ugent.psb.setas.bdmodel.oldclasses;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;

import be.ugent.psb.setas.bdmodel.model.FindMaxArray;
import be.ugent.psb.setas.bdmodel.model.LikeLihood;
import be.ugent.psb.setas.bdmodel.model.MathOperations;
import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.bdmodel.parsers.NewickParser;
import be.ugent.psb.setas.bdmodel.parsers.WGMparser;

public class Interpolation {

	public Node root;
	public ArrayList<Node> nodes;
	public int rootSize;
	int counter;

	public double optimizeLambda(double minLam, double maxLam, double precision) {

		double optimalLambda;
	
		FindMaxArray fm = new FindMaxArray();

		double rangeLenx = maxLam - minLam;
		
		int numberOfPartitions = (int) (rangeLenx / precision);
		
		double[] testLambdas = new double[numberOfPartitions + 1];
			
		double[] correspondingLk = new double[testLambdas.length];
		double[] correspondingLogLk = new double[testLambdas.length];
		
		for (int i = 0; i < numberOfPartitions + 1; i++) {
			testLambdas[i] = minLam + (i * precision);
			
			/* to include angiosperms */
			LikeLihood lk = new LikeLihood(testLambdas[i],root.getmaxNodeSize());
			double [] crrsLks= lk.calcInternalLk(nodes);
			correspondingLk[i] = crrsLks[rootSize];			
//			correspondingLk[i] = lk.calcLkForRootSize(root, rootSize, testLambdas[i]);

		}
				
		correspondingLogLk = MathOperations.giveLogArray(correspondingLk);

		int maxIndex = fm.findIndexOfMaxValue(correspondingLogLk);	
		
		optimalLambda = testLambdas[maxIndex];
		
		precision = precision * 0.1;
		
		counter +=1;
		
		if (precision > 0.001 && numberOfPartitions >=1) {

			double newMaxLam = optimalLambda + (double) (1.00/counter);
			double newMinLam = optimalLambda - (double)(1.00/counter);
			return(optimizeLambda(newMinLam, newMaxLam, precision));
		}
		
		else{
		return optimalLambda;}
	}

//	public static void main(String[] args) {
//
//		long startTime = System.currentTimeMillis();
//
//		Log lg = new Log();
//		FindMaxArray fma = new FindMaxArray();
//		NewickParser np = new NewickParser();
//		Node root = np
//				.buildTree("/home/setas/workspace/BirthDeathModel/src/Files/EudicotTree-Arab");
//
//		WGMparser wgm = new WGMparser();
//		List<List<String>> wgdList = wgm
//				.readInputFile("/home/setas/workspace/BirthDeathModel/src/Files/WGDeudicotTree-Arab");
//		root.insertWGM(root, wgdList);
//
//		ArrayList<Node> leaves = root.getLeaves();
//		root.setLeaves(leaves);
//		root.setNumberOfLeaves(leaves.size());
//
//		int[] tpp = new int[root.getNumberOfLeaves()];
//
//		tpp[0] = 7;
//		tpp[1] = 22;
//		tpp[2] = 12;
//		tpp[3] = 7;
//		tpp[4] = 10;
//		tpp[5] = 7;
//		tpp[6] = 11;
//		tpp[7] = 5;
//		tpp[8] = 5;
//		tpp[9] = 9;
//
//		root.setLeafValues(tpp);
//
//		LikeLihood lk = new LikeLihood();
//
//		Queue<Node> queue = lk.PostOrder(root);
//		ArrayList<Node> aryln = lk.qToArrayNode(queue);
//
//		Interpolation inter = new Interpolation();
//		inter.root=root;
//		inter.nodes = aryln;
//		inter.rootSize = 5;
//
//		double bestLambda = inter.optimizeLambda(0, 10.3, 1);
//		
//		long endTime = System.currentTimeMillis();
//		System.out.println("Took " + (endTime - startTime) + " mili seconds");
//
//	}

}
