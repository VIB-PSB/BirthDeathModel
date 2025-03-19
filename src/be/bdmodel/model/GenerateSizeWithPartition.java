package be.ugent.psb.setas.bdmodel.model;

//
//import java.util.Map;
//import java.util.TreeMap;
//
//import be.ugent.psb.setas.bdmodel.model.ProbCalculator;
//
//
public class GenerateSizeWithPartition {
//
//
//	static double deltaT = 0.1;
//
//	/**
//	 * Given the size of a Parent node, generates a size for its child node,
//	 * sampled from a probablity distribution described in the class
//	 * Probcalculator, via MCMC method.
//	 * 
//	 * Parent nodes can not take gene count = 0 since the be.ugent.psb.setas.bdmodel.model does not support
//	 * 0 values in the internal nodes: absorbing state of the system
//	 *
//	 * @param maxNodeSize
//	 * @param parentSize gene counts at the parent node
//	 * @param lambda gene duplication rate
//	 * @param wholeBrLenx branch length
//	 * @return size of the child node
//	 */
//	
//	
//	public static int generateSize_partitioning(int maxNodeSize,
//			int parentSize, double lambda, double wholeBrLenx) {
//		
////		System.out.println("Warning: generateSize_partitioning working with no cache");
//	
//		return generateSize_partitioning(maxNodeSize,parentSize, lambda, wholeBrLenx, null);
//	}
//	
//	public static int generateSize_partitioning(int maxNodeSize,
//			int parentSize, double lambda, double wholeBrLenx, ProbCalculator probCache) {
//
//
//		GenerateSize gs = new GenerateSize(probCache);
//		
//		if (wholeBrLenx > deltaT) {
//			double counter = 1;
//			int numOfPartitions = (int) (Math.floor(wholeBrLenx / deltaT));
//			double remainingBlen = wholeBrLenx - (numOfPartitions) * deltaT;
//			
//			int childSizeOld = gs.generateSize(maxNodeSize,parentSize, lambda, deltaT);
//			int childSizeNew = childSizeOld;
//
//			while (counter <= numOfPartitions) {
//
//				childSizeNew = gs.generateSize(maxNodeSize,childSizeOld, lambda, deltaT);
//				
//				childSizeOld = childSizeNew;
//				counter += 1;
//			}
//
//			if (Math.abs(remainingBlen - 0) < 1e-10) {
//				return childSizeNew;
//			}
//
//			else {
//				int childSizeLast = gs.generateSize(maxNodeSize,
//						childSizeNew, lambda, remainingBlen);
//				return childSizeLast;
//			}
//		}
//
//		else {
//			return (gs.generateSize(maxNodeSize, parentSize, lambda,wholeBrLenx));
//		}
//	}
//
//	/*
//	 * Because at the leaves, the gene family can go extinct, so we can generate
//	 * 0 for leaf sizes : j=0
//	 */
//	public static int generateSizeForleaves_partitioning(int maxNodeSize,
//			int parentSize, double lambda, double wholeBrLenx) {
//		
//		System.out.println("Warning: generateSizeForLeaves_partitioning working with no cache");
//		return generateSizeForleaves_partitioning(maxNodeSize, parentSize,lambda, wholeBrLenx, null);
//	}
//
//	public static int generateSizeForleaves_partitioning(int maxNodeSize,
//			int parentSize, double lambda, double wholeBrLenx, ProbCalculator probCache) {
//		
//		GenerateSize gs = new GenerateSize(probCache);
//		if (wholeBrLenx > deltaT) {
//		double counter = 1;
//
//		int numOfPartitions = (int) (Math.floor(wholeBrLenx / deltaT));
//
//		double remainingBlen = wholeBrLenx - (numOfPartitions) * deltaT;
//
//		int childSizeOld = gs.generateSize(maxNodeSize, parentSize,
//				lambda, deltaT);
//		
//		int childSizeNew = childSizeOld;
//
//		// only run for one to the end interval
//		while (counter < numOfPartitions) {
//
//			childSizeNew = gs.generateSize(maxNodeSize, childSizeOld,
//					lambda, deltaT);
//			childSizeOld = childSizeNew;
//			counter += 1;
//		}
//
//		// you are already at the leaf:
//		if (Math.abs(remainingBlen - 0) < 1e-10) {
//
//			childSizeNew = gs.generateSizeForleaves(maxNodeSize,
//					childSizeOld, lambda, deltaT);
//			return childSizeNew;
//		}
//
//		else { // there is still one more step left to reach the leaf
//
//			childSizeNew = gs.generateSize(maxNodeSize, childSizeNew,
//					lambda, remainingBlen);
//			int childSizeLast = gs.generateSizeForleaves(maxNodeSize,
//					childSizeNew, lambda, remainingBlen);
//
//			
//			return childSizeLast;
//		}
//		}
//		else {
//			int childSizeLast = gs.generateSizeForleaves(maxNodeSize, parentSize, lambda,
//					wholeBrLenx);	
//			return childSizeLast;
//		}	
//	}
//
////	public static void main(String[] args) {
////
////		// for(int i=0; i<10000; i++){
////		int d = GenerateSizeWithPartition.generateSize_partitioning(100, 1, 10,
////				0.831625);
////
////		// if(d < 1 || d > 100){
////		// System.out.println(d);
////		// }
////`+
////		// System.out.println("d "+ d);
////		// }
////
////	}
//
}
