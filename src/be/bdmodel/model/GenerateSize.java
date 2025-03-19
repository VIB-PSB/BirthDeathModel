package be.ugent.psb.setas.bdmodel.model;

import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import be.ugent.psb.setas.bdmodel.model.ProbCalculator;

public class GenerateSize {

	/**
	 * given the size of a Parent node, generates a size for its child node, sampled
	 * from a probablity distribution described in the class Probcalculator, via
	 * MCMC method
	 * 
	 * @param parentSize = gene counts at the parent node
	 * @param lambda     = gene duplication/loss rate
	 * @param t          = branch length
	 * @return
	 */

	private ProbCalculator probCalc;
	private int lengthOfMCMC;
	private int minNodeSize;
	private int maxNodeSize;

	public GenerateSize() {
		this.probCalc = new ProbCalculator();
		this.lengthOfMCMC = 5000;
	}

	public GenerateSize(ProbCalculator probCache) {
		this.probCalc = probCache;
		this.lengthOfMCMC = 5000;

	}

	public GenerateSize(ProbCalculator probCache, int lMCMC) {
		this.probCalc = probCache;
		this.lengthOfMCMC = lMCMC;
	}
	
	public GenerateSize(ProbCalculator probCache, int lMCMC, int minNodeSize,int maxNodeSize) {
		this.probCalc = probCache;
		this.lengthOfMCMC = lMCMC;
		this.minNodeSize=minNodeSize;
		this.maxNodeSize = maxNodeSize; 
	}
	
	/*
	 * Parent nodes can not take gene count =0 since the
	 * be.ugent.psb.setas.bdmodel.model does not support 0 values in the internal
	 * nodes: absorbing state of the system
	 */

	public int generateSize(int parentSize, double lambda, double t)  {
	
        parentSize= bringInRange(minNodeSize, maxNodeSize, parentSize);
		Random random = new Random();
		
		/* So if the parentSize does not change form one step to the next we don't re-calculate it*/
		double probNoChange = this.probCalc.probCalc(lambda, t, parentSize, parentSize);
		int nextGeneratedSize = parentSize;

		/* 1000 iterations of MCMC */
		for (int k = 1; k <= lengthOfMCMC ; k++) {

			 /*random = from + rndGenerator.nextInt(to - from + 1) */       
		     int proposedChildSize = 1+random.nextInt(maxNodeSize) ;

		     double newProb=this.probCalc.probCalc(lambda, t, nextGeneratedSize, proposedChildSize);
		     if(newProb > probNoChange) {nextGeneratedSize = proposedChildSize;}
		     else {
		    	 
		     double fraction = newProb / probNoChange;
			
			 double alpha = Math.min(1.0, fraction);
//			 double randomValue = rangeMin + (rangeMax - rangeMin) * r.nextDouble();
			 double u = random.nextDouble();

			 if (alpha >= u) {
				nextGeneratedSize = proposedChildSize;
				probNoChange = this.probCalc.probCalc(lambda, t, nextGeneratedSize,nextGeneratedSize );
			 }		
		    }
//			else{Do nothing, keep same parent size, generate new random child size and go on with that..}
		}
		return nextGeneratedSize;
	}

	/*
	 * Because at the leaves, the gene family can go extinct, so we can generate 0
	 * for leaf sizes : c=0
	 */
	public int generateSizeForleaves(int parentSize, double lambda, double t) {

		Random random = new Random();

        parentSize= bringInRange(minNodeSize, maxNodeSize, parentSize);

		double probNoChange = this.probCalc.probCalc(lambda, t, parentSize, parentSize);
		int nextGeneratedSize = parentSize;
		double u = random.nextDouble();

		for (int k = 1; k < lengthOfMCMC + 1; k++) {

			int proposedChildSize = 1+random.nextInt(maxNodeSize); // because lk array is hard wired for 100 slots

			double newProb = this.probCalc.probCalc(lambda, t, nextGeneratedSize, proposedChildSize);

			if (newProb > probNoChange) {
				nextGeneratedSize = proposedChildSize;
			} else {
				double fraction = newProb / probNoChange;

				double alpha = Math.min(1.0, fraction);

				if (alpha >= u) {

					nextGeneratedSize = proposedChildSize;

//					/*
//					 * If you produced 0 counts for some leaf you can not go on in the MCMC chain
//					 * and the chain will stop and you return 0
//					 */
//					if (proposedChildSize == 0) {
//						return 0;
//					}
//
//					else {
//						probNoChange = this.probCalc.probCalc(lambda, t, nextGeneratedSize, nextGeneratedSize);
//					}
				}
			}
		}

		return nextGeneratedSize;

	}
	
	public int bringInRange( int minNodeSize, int maxNodeSize, int prob) {
		
		if(prob < minNodeSize) {return minNodeSize;}
		else if(prob > maxNodeSize) {return maxNodeSize;}
		else {
			return prob;
		}
		
	}

//	public static void main(String[] args) {
//
//		GenerateSize gs = new GenerateSize();
//		
//		for(int i=0; i<10000; i++){
//			
//		int d = gs.generateSize(100, 1, 1.0, 0.8);
//		
//		if(d < 1 || d > 100){
//			System.out.println("error: "+ d);
//		}
//		}
//		
//	}

}
