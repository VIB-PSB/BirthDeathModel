package utils.bdmodel;

import java.util.Random;

public class GenerateSize {

	/**
	 * given the size of a Parent node, generates a size for its child node, sampled
	 * from a probability distribution described in the class ProbCalculator, via
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
		this.minNodeSize = minNodeSize;
		this.maxNodeSize = maxNodeSize; 
	}

        /*
	 * Parent nodes can not take gene count =0 since the
	 * be.ugent.psb.setas.bdmodel.model does not support 0 values in the internal
	 * nodes: absorbing state of the system
	 */

        
	public int generateSize(int parentSize, double lambda, double t)  {
            
            Random random = new Random();

            parentSize= bringInRange(minNodeSize, maxNodeSize, parentSize);

            double probNoChange = this.probCalc.probCalc(lambda, t, parentSize, parentSize);
            int currentSize = parentSize;
	
       
            for (int k = 0; k < lengthOfMCMC; k++) {

                int nextstep;
        
                double v = random.nextDouble();
                if(v<0.5){
                    nextstep = -1;
                }
                else{
                    nextstep = 1;
                }
            
                // No internal node sizes < 1 and no internal node sizes > maxNodeSize
                if ((currentSize == 1 && nextstep == -1) || (currentSize == maxNodeSize && nextstep == 1)){
                    continue;
                }
                else{
                    int proposedChildSize = currentSize + nextstep;
                    double newProb = this.probCalc.probCalc(lambda, t, parentSize, proposedChildSize);
                    double alpha = newProb / probNoChange;
                    double u = random.nextDouble();
                    if (u <= alpha) {
                        currentSize = proposedChildSize;
                        probNoChange = this.probCalc.probCalc(lambda, t, parentSize, currentSize);
                    }
                    // else step is not accepted and currentSize remains the same                
                }
            }
        
            return currentSize;
        }
        
	public int generateSizeForleaves(int parentSize, double lambda, double t) {

            Random random = new Random();

            parentSize= bringInRange(minNodeSize, maxNodeSize, parentSize);

            double probNoChange = this.probCalc.probCalc(lambda, t, parentSize, parentSize);
            int currentSize = parentSize;
	
       
            for (int k = 0; k < lengthOfMCMC; k++) {

                int nextstep;
        
                double v = random.nextDouble();
                if(v<0.5){
                    nextstep = -1;
                }
                else{
                    nextstep = 1;
                }
            
                // Leaf sizes = 0 are problematic as this.probCalc.probCalc(lambda, t, currentSize, currentSize) doesn't work for transition of 0 to 0 genes
                // Assume no zero counts at leaves (or biologically, 'conserved' gene families do not get lost completely)
                // No leaf sizes < 1 and no leaf sizes > maxNodeSize
                
                if ((currentSize == 1 && nextstep == -1) || (currentSize == maxNodeSize && nextstep == 1)){
                    continue;
                }
                else{
                    int proposedChildSize = currentSize + nextstep;
                    double newProb = this.probCalc.probCalc(lambda, t, parentSize, proposedChildSize);
                    double alpha = newProb / probNoChange;
                    double u = random.nextDouble();
                    if (u <= alpha) {
                        currentSize = proposedChildSize;
                        probNoChange = this.probCalc.probCalc(lambda, t, parentSize, currentSize);
                    }
                    // else step is not accepted and currentSize remains the same                
                }
            }
        
            return currentSize;
        }
            
                    
	
	public int bringInRange( int minNodeSize, int maxNodeSize, int prob) {
		
		if(prob < minNodeSize) {return minNodeSize;}
		else if(prob > maxNodeSize) {return maxNodeSize;}
		else {
			return prob;
		}
		
	}
}
