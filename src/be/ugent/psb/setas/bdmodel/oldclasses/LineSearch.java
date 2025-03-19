package be.ugent.psb.setas.bdmodel.oldclasses;

import java.util.ArrayList;

import be.ugent.psb.setas.bdmodel.model.FindMaxArray;
import be.ugent.psb.setas.bdmodel.model.LikeLihood;
import be.ugent.psb.setas.bdmodel.model.MathOperations;
import be.ugent.psb.setas.bdmodel.model.Node;


public class LineSearch {

	/**
	 * @param args
	 */
	
	public Node root;
	public int numOfPartitions;
	
	public double lineSearch(Node root, int rootSize){
		
		double [] lambdas = new double [numOfPartitions];
		double [] likelihoods = new double [numOfPartitions];
		double [] logLikelihoods=  new double[numOfPartitions]; 
		ArrayList<Node> arln = root.postOrder(root);
		
		
		double increment= 0.0001;
		
		for(int i=1; i<numOfPartitions; i++){
			
			lambdas[i] = i* increment;
			
			LikeLihood lk=new LikeLihood(lambdas[i],root.getmaxNodeSize()+1);
			double [] lks = lk.calcInternalLk(arln);
			likelihoods[i]= lks[rootSize]; 
				
//			System.out.println("likelihood for lambda= " +lambdas[i] +" = "+ likelihoods[i]);
			
		}
		
		//define another vector that stores the maximum likelihood for each row, i.e for each root size
		// also print the related lambda
		
		logLikelihoods = MathOperations.giveLogArray(likelihoods) ;
		FindMaxArray fma = new FindMaxArray(); 
		
		double maxLk= fma.findMaxValue(likelihoods);
//		System.out.println("maximum likelihood: "+ maxLk);
		
		double maxLogLk= fma.findMaxValue(logLikelihoods);
//		System.out.println("maximum log-likelihood: "+ maxLogLk);
		
		int maxIndex = fma.findIndexOfMaxValue(logLikelihoods);
		
//		System.out.println("max Index: "+ maxIndex);
		
		double bestLam= lambdas[maxIndex];
		
		return bestLam;
	}
	
//	public static void main(String[] args) {
//		// TODO Auto-generated method stub
//
//	}

}
