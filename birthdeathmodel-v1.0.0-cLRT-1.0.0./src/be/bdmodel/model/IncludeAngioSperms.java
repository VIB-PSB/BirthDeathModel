package be.ugent.psb.setas.bdmodel.model;
//
public class IncludeAngioSperms {
//
//	/*
//	 * 
//	 * To include the angiosperms and the WGT in the ancestor of Eudicots, we can just
//	 * change the likelihood array accordingly. 
//	 */
//	
//	public static double[] convertLkArray (double [] oldLks, double lambda, double wgtBlen, double oldrootBlen, int factor){
//		return IncludeAngioSperms.convertLkArray(oldLks, lambda, wgtBlen, oldrootBlen, null,factor);
//	}
//	
//	public static double[] convertLkArray (double [] oldLks, double lambda, double wgtBlen, double oldrootBlen, ProbCalculator probCache, int factor){
//		
//		int numOfTestRootSizes = oldLks.length;
//		
//		double newEudicotBlen = oldrootBlen - wgtBlen;
//		
//		double [] afterWGTLks = new double [numOfTestRootSizes];
//		
//		ProbCalculator probCalc = probCache;
//		
//		for(int j=factor; j< numOfTestRootSizes; j=j+factor){
//			
//			double r=0;
//			
//			for (int k = 1; k < numOfTestRootSizes; k++) {
//				r += probCalc.probCalc(lambda, newEudicotBlen, j, k)* oldLks[k];
//			}
//			
//			afterWGTLks[j] = r; 
//		}
//		
//		double [] angioLiks = new double[numOfTestRootSizes];
//				
//		for (int j = 1; j < numOfTestRootSizes; j++) {
//
//			double r = 0;
//
//			for (int k = factor; k < numOfTestRootSizes; k = k+ factor) {
//				r += probCalc.probCalc(lambda, wgtBlen, j, k / factor)
//						* afterWGTLks[k];
//			}
//
//			angioLiks[j] = r ;
//		}
//			
//		return angioLiks;
//	}
//	
//	public static int generateSizeAtEudicots(int maxNodeSize, int s, double lambda, double wgtBlen, double oldrootBlen){
//		
////		System.out.println("Warning: IncludeAngioSperms is using no Cache");
//		return generateSizeAtEudicots(maxNodeSize, s, lambda, wgtBlen, oldrootBlen,null);
//	}
//	
//	public static int generateSizeAtEudicots(int maxNodeSize, int s, double lambda, double wgtBlen, double oldrootBlen, ProbCalculator probCache){
//		
//		int eudicotsSize = 0;
//		
//		/* For all the other nodes, because they are included in the structure of the tree, we use partitionBranches. Hence all branches are already chopped up into
//		 * piences of lenght deltaT. So we can use the normal generateSize method to generate a random siz for the child node. 
//		 * However since the last branch, connecting Eudicots to Angiosperms is not added to the structure of the tree, we need to be careful to apply the same 
//		 * chopping up step before we generate size. For instance if we use lambda_max =10, then we have to make sure to use GenSizeWithPartitioning to generate a size, 
//		 * recursively for pieces of the branch having size 0.1 So we use GenSizePartitioning which was initially only used for SC project.
//		 * 
//		 * It is not important to integrate these new changes to the tree, i.e. to add VirtualNodes on this branch, since the structure is simple and unique and we are
//		 * treating it individually anyway.*/
//		
//		int beforeWGT = GenerateSizeWithPartition.generateSize_partitioning(maxNodeSize, s, lambda, wgtBlen, probCache);
//		int afterWGT = 3* beforeWGT;
////		int afterWGT = 2*beforeWGT; // for the tree of monocots only 
//		
//		double newEudicotBlen = oldrootBlen - wgtBlen;
//		
//		eudicotsSize = GenerateSizeWithPartition.generateSize_partitioning(maxNodeSize, afterWGT, lambda, newEudicotBlen, probCache); 
//		
//		return eudicotsSize;
//		
//		
//	}
}
