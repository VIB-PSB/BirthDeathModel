package utils.bdmodel;

import java.util.ArrayList;

// TODO: remove un-necessary import cause now in same bdmodel directory?
// import utils.bdmodel.CuttingPlaneMethod;
// import utils.bdmodel.Node;
// import utils.bdmodel.ProbCalculator;

public class EstimateLambdaRootSize {
	
	final double stepSize = 1e-4; // Step calculating derivative
	final double deltaLocalMoves = 1e-1;
	final double tolD = 1e-3;
	final double tolF = 1e-4;
	final double minInterval = 1e-2;
	final double maxInterval = 9.99;
	final double precisionLambda = 1e-5;
	
	// TODO: why says it's not used while in original code yes?
	private int maxRootNodeSize;
	
	public EstimateLambdaRootSize(int rootSize) {
		this.maxRootNodeSize = rootSize;
	}
	
	public double[] cpm_optimize_r_lambda(ArrayList<Node> speTree, Node root, ProbCalculator probCalc) {
		
		double[] rStar_lambdaStar_LoglkStar = new double [3];
		
		for(int testRootSizes = 1; testRootSizes <= 10; testRootSizes++) {
			
			root.setValue(testRootSizes);
			
			double[] tmp= cpm_fixedRootSize (speTree, testRootSizes, root, probCalc);
			
			if(testRootSizes==1) {
				rStar_lambdaStar_LoglkStar = tmp;}
			
			else { // Find the rootSize with maximum loglk
				if(tmp[2] > rStar_lambdaStar_LoglkStar[2]) {
					rStar_lambdaStar_LoglkStar = tmp;			
				}
			}
		}
		
		return rStar_lambdaStar_LoglkStar;
	}
	
	public double[] cpm_fixedRootSize (ArrayList<Node> speTree, int rootSize, Node root, ProbCalculator probCalc){
		
		double[] rStar_lamStar_LoglkStar = new double [3];
		rStar_lamStar_LoglkStar[0] = rootSize; 
		
		CuttingPlaneMethod cpm = new CuttingPlaneMethod(speTree,
		rootSize, stepSize, deltaLocalMoves, tolD, tolF,
		minInterval, maxInterval,
		precisionLambda, probCalc, root);	 
		cpm.findOptimalLambda();	
		rStar_lamStar_LoglkStar[1] = cpm.getOptimalLambda();
		rStar_lamStar_LoglkStar[2] = cpm.getFoptimalLambda();
		return rStar_lamStar_LoglkStar;
	}
	
}
