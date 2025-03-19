package be.ugent.psb.setas.bdmodel.model;

import java.util.ArrayList;
//import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import be.ugent.psb.setas.bdmodel.parsers.NewickParser;
import be.ugent.psb.setas.bdmodel.parsers.ReadGFcountsFile;
import be.ugent.psb.setas.bdmodel.parsers.WGMparser;
//import be.ugent.psb.setas.bdmodel.parsers.BuildTree;

public class Main {
	
	public static void main(String[] args) {
		
		long startTime = System.currentTimeMillis();
			
		/* The gene counts file should have the leaves in the order of the tree, or change the order in the tree*/
		ReadGFcountsFile rgf = new ReadGFcountsFile();
		List<List<Integer>> tbl = rgf.read_unique(args[2]);
		ArrayList<String> gfIDs = rgf.getGfIDs_unique();

		List<List<Integer>> counts = tbl;		
		
		/* which column = which gf counts */
		int gf = Integer.parseInt(args[3]);	
				
		List<Integer> li = counts.get(gf);			

		int[] originalObservation = new int[li.size()];

		for (int m = 0; m < li.size(); m++) {
			originalObservation[m] = li.get(m);
		}
		
		int maxGeneCount = FindMaxArray.findMaxValueIntArray(originalObservation); 
		int maxNodeSize = Math.max(100,maxGeneCount);
			
		NewickParser np = new NewickParser();
		Node root = np.buildTree(args[0],maxNodeSize);

		WGMparser wgm = new WGMparser();
		List<List<String>> wgdList = wgm.readInputFile(args[1]);
		root.insertWGMsToTheTree(wgdList);

		ArrayList<Node> leaves = root.getLeaves();
//		root.setLeaves(leaves);
//		root.setNumberOfLeaves(leaves.size());

		root.setLeafValues(originalObservation);
		
	    PartitionBranches pb = new PartitionBranches(root,0.1);
		pb.insertAllVNSonAllBranches();
		
		Queue<Node> queue2 = root.postOrder();
		ArrayList<Node> aryln = new ArrayList<Node>(queue2);
		
		for(int i=0; i<aryln.size();i++){
			System.out.println(aryln.get(i).getName());
		}

		double supremumLambda = 10.0;		
		double initialStepSize = 0.1;

		/* number of rows = number of lambdas */
		int rows = (int) Math.floor((supremumLambda / initialStepSize));

		/* number of columns = number of root sizes */
		int columns = 21;
		
		Matrix mt = new Matrix();
		double[][] d = mt.calcMatrix(aryln, rows, columns, initialStepSize);

		double[] maxLambdas = new double[columns];
		double[] maxLogLk = new double[columns];
//		double[] pValues = new double[columns];

//		Pvalues pv = new Pvalues(root, aryln, 1000);
		/*
		 * The root Node is Eudicot, 
		 * We changed the likelihood array to calculate lks at the angiosperms.
		 */
	
//		double wgtBlen = 0.260;
//		double oldrootBlen = 0.275776;

		/* means size at angiosperms is s */
		for (int s = 1; s < columns; s++) {

			/* take one column of matrix corresponding to one root size */
			double[] lkForOneS = new double[rows];
			double[] lgLkForOneS = new double[rows];

			/* get rid of rows or lambdas: */
			for (int lam = 1; lam < rows; lam++) {
				lkForOneS[lam] = d[lam][s];
			}

			// we are already taking into account includeAngioSperm inside Likelihood class, the following eudicots size story is only for P-val calcs.
			lgLkForOneS = MathOperations.giveLogArray(lkForOneS);
			maxLogLk[s] = FindMaxArray.findMaxValue(lgLkForOneS);

			int maxIndex = FindMaxArray.findIndexOfMaxValue(lgLkForOneS);
			maxLambdas[s] = maxIndex * initialStepSize;

//			int eudicotSize = IncludeAngioSperms.generateSizeAtEudicots(root.getmaxNodeSize(), s, maxLambdas[s],
//					wgtBlen, oldrootBlen);
//
//			/*If we are using the whole tree, rooted at angiosperms.*/
////			int eudicotSize = s;
//			pValues[s] = pv.calculateConditionalPvalues(s, eudicotSize,
//					maxLambdas[s], maxLogLk[s]);
			
		}

		
		double maxLogLkOverAllS = FindMaxArray.findMaxValue(maxLogLk);
		
		int bestRootSize = FindMaxArray.findIndexOfMaxValue(maxLogLk);	

//		double supPvalue = fma.ValueFindMax(pValues);
//			
		System.out.print(gfIDs.get(gf)+"\t"+bestRootSize + "\t" + maxLambdas[bestRootSize]
		+ "\t" + maxLogLkOverAllS);
//		System.out.print(gfIDs.get(gf)+"\t"+bestRootSize + "\t" + maxLambdas[bestRootSize]
//				+ "\t" + maxLogLkOverAllS + "\t" + pValues[bestRootSize]+"\t"
//				+ supPvalue);
//		
//		 if (supPvalue > 0.01) {
//		 System.out.print("\t"+"follows");
//		
//		 } else {
//		 System.out.print("\t"+"NOT-follows");
//		 }
//		 
//		for(int j=1; j<columns;j++){
//			
//			System.out.print("\t"+maxLambdas[j]);
//			System.out.print("\t" + pValues[j]);
//			
//		}
//		
		long endTime = System.currentTimeMillis();
		System.out.print("\t" + (endTime-startTime));

		/******************************************/
		// EM em = new EM();
		// em.epsilon = 0.001;
		// em.root = root;
		// double bestLambda = em.emForLambda(4);
		// int bestRootSize = root.getOptimalSize();
		// double[] bestLik = lk.calcInternalLk(aryln, bestLambda);
		// double [] lgBestLik = lg.giveLogArray(bestLik);
		// System.out.println("Final: Root Size: "+
		// bestRootSize+" Lambda: "+bestLambda+" LogLk: "+
		// lgBestLik[bestRootSize]);

		/***************************************/



	}	
}
