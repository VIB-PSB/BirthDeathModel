package be.ugent.psb.setas.bdmodel.model;

import java.util.ArrayList;

public class Matrix {

	/**
	 * calculates likelihoods for all lambdas between 0 and 1, increasing by
	 * stepSize every time.(rows) and for all root sizes (= columns)
	 * 
	 * @param args
	 */

	public double[][] calcMatrix(ArrayList<Node> arln, int rows, int columns,
			double stepSize) {

		int numOfTestRootSizes = arln.get(0).getmaxNodeSize() + 1;

		/*
		 * int columns = number of root sizes; int rows = number of partitions
		 */
		double[][] m = new double[rows][columns];
		double[] temp = new double[columns];

		for (int lam = 0; lam < rows; lam++) {

			double lambda = lam * stepSize;
			LikeLihood lk = new LikeLihood(lambda, numOfTestRootSizes);
			temp = lk.calcInternalLk(arln);

			for (int s = 1; s < columns; s++) {

				m[lam][s] = temp[s];

			}
		}
		return m;
	}

//	 public static void main(String[] args) {
//	
//		 
//		 Matrix mt = new Matrix();
//		 
//		 double [][] m = mt.calcMatrix(arln, 3000, 3, 0.01);
//		 
//		 
//	
//	 }

}
