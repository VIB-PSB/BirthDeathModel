package be.ugent.psb.setas.bdmodel.oldclasses;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import be.ugent.psb.setas.bdmodel.parsers.NewickParser;
import be.ugent.psb.setas.bdmodel.parsers.ReadGFcountsFile;
import be.ugent.psb.setas.bdmodel.parsers.WGMparser;

public class EucDistance {

	public double calcEucDistance(int[] a, int[] b) {

		double sum = 0;
		double dist = 0;

		if (a.length == b.length) {

			for (int i = 0; i < a.length; i++) {

				sum += Math.pow(a[i] - b[i], 2);
			}

			dist = Math.sqrt(sum);
		}

		else {
			System.out
					.println("Euclidean-dist can't be calculated: length of arrays are different");
		}
		return dist;

	}

//	public static void main(String args[]) {
//
//		NewickParser np = new NewickParser();
//		Node root = np
//				.buildTree("/home/setas/workspace/BirthDeathModel/src/Files/EudicotTree");
//
//		WGMparser wgm = new WGMparser();
//		List<List<String>> wgdList = wgm
//				.readInputFile("/home/setas/workspace/BirthDeathModel/src/Files/EudicotTreeWGM");
//		root.insertWGMsToTheTree(wgdList);
//
//		ArrayList<Node> leaves = root.getLeaves();
//		root.setLeaves(leaves);
//		root.setNumberOfLeaves(leaves.size());
//
//		GenerateObservations gobs = new GenerateObservations();
//
//		int[] idealObservations = gobs.generateObservation(root, 2, 0);
//		int[] includeWGTidealObs = new int[idealObservations.length];
//
//		for (int k = 0; k < includeWGTidealObs.length; k++) {
//
//			includeWGTidealObs[k] = 3 * idealObservations[k];
//		}
//
//		ReadGFcountsFile rgf = new ReadGFcountsFile();
//		List<List<Integer>> tbl = rgf
//				.readInputFile("/home/setas/workspace/BirthDeathModel/src/Files/GF_Eudicots");
//
//		List<List<Integer>> counts = rgf.removeColumnsWithZero(tbl);
//		int numberOfGfs = counts.size();
//		
//		 int[] tpp = new int[root.getNumberOfLeaves()];
//		 tpp[0] = 7;
//		 tpp[1] = 22;
//		 tpp[2] = 12;
//		 tpp[3] = 7;
//		 tpp[4] = 10;
//		 tpp[5] = 7;
//		 tpp[6] = 11;
//		 tpp[7] = 10;
//		 tpp[8] = 10;
//		 tpp[9] = 5;
//		 tpp[10] = 5;
//		 tpp[11] = 9;
//
//		/* which column = which gf counts */
//		// int gf = Integer.parseInt(args[3]);
//
//		EucDistance eud = new EucDistance();
//		double[] distances = new double[numberOfGfs];
//		double[] temp = new  double [numberOfGfs];
//		
//	
////		int gf= 45;
////		List<Integer> li = counts.get(gf);
////
////		int[] originalObservation = new int[li.size()];
////
////		for (int m = 0; m < li.size(); m++) {
////		originalObservation[m] = li.get(m);}
//		
//		
//		distances[0] = eud.calcEucDistance(includeWGTidealObs,
//				tpp);
//		
////		System.out.println(distances[0]);
//		
//
////
//////			distances[gf] = eud.calcEucDistance(includeWGTidealObs,
//////					originalObservation);
////			System.out.println(gf+"  "+distances[gf]);
////			
////			//make a copy of this array to sort later:
////			temp[gf] = distances[gf];
////		}	
//
////		/* to sort the results with indexes: */
////		
////		Arrays.sort(temp); // sort ascending
////		
//////		for(int i=0; i<temp.length; i++){
//////			
//////			System.out.println("temp "+temp[i]);
//////		}
////					
////		// final array of indexes
////		int index_array[] = new int[numberOfGfs];
////
//////		boolean found = false;
////		
//////		boolean [] visited = new boolean [numberOfGfs];
////		
////		 /*iteretate on the sorted array */
////		for (int i = 0; i < numberOfGfs; i++) {
////			
//////			System.out.println(temp[i]);
////			
////			//iterate on the original array
////			for(int j=0; j< numberOfGfs; j++){
////							
////				if(Math.abs(distances[j] - temp[i]) == 0){
////					
////					System.out.println("j: "+ j);
//////					visited[j]= true;
////					index_array[i]= j;
////				}
////			
////				
////			}			
////
////		}
////
////		for (int i = 0; i < numberOfGfs; i++) {
////
////			System.out.println(index_array[i]+"\t"+distances[index_array[i]]);
////	
////		}
//
//	}

}
