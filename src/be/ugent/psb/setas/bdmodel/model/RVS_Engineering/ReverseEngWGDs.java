package be.ugent.psb.setas.bdmodel.model.RVS_Engineering;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.Scanner;

import be.ugent.psb.setas.bdmodel.model.Branch;
import be.ugent.psb.setas.bdmodel.model.ChooseKoutOfN;
import be.ugent.psb.setas.bdmodel.model.FindMaxArray;
import be.ugent.psb.setas.bdmodel.model.LikeLihood;
import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.bdmodel.model.ProbCalculator;
import be.ugent.psb.setas.bdmodel.parsers.NewickParser;
import be.ugent.psb.setas.bdmodel.parsers.ReadGFcountsFile;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;
import be.ugent.psb.setas.bdmodel.parsers.WGMparser;

public class ReverseEngWGDs {

	private Node rootOfTree;
	private ArrayList<Node> arln;
	
	public Node parentToWGD;
	public Node childToWGD;
	/* numPartitionsOnSusBr must be at least 2, because otherwise th whole branch is treated equally, not knowing where to put WGM, how to define its brlen, etc
	 * If numPartitions=2 then there's only one place: in the middle. To be.ugent.psb.setas.bdmodel.test more possibilities, partition it to 3, 4, etc*/
	public int numPartitionsOnSusBr;

	private String mgfFile;
	private String treeFile;
	private ProbCalculator probCache;
	private List<List<String>> fixedWGDList;

	public ReverseEngWGDs(String treeFile, String mgfFile, Node rootOfTree,
			ProbCalculator probCache) {
		this.arln = rootOfTree.postOrder(rootOfTree);
		this.treeFile = treeFile;
		this.mgfFile = mgfFile;
		this.rootOfTree = rootOfTree;
		this.probCache = probCache;
	}

	/**
	 * The MGF file should be a file in each line one Marker-Gene-Family ID and
	 * tab separated: rootSize*, lambda*, and then the gene counts
	 */
	public List<List<Double>> readMGFdata() {

		FileReader fin = null;
		try {
			fin = new FileReader(mgfFile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		/* The first line is a header: */
		Scanner sc = new Scanner(fin);
		sc.nextLine();

		List<List<Double>> MGFtable = new ArrayList<List<Double>>();

		while (sc.hasNextLine()) {

			String line = sc.nextLine();

			String[] numbers = line.split("\t");

			/** Now we have the numbers as strings in the array "numbers", so: */
			List<Double> nums = new ArrayList<Double>();

			for (int i = 1; i < numbers.length; i++) { // i=0 ~ GF-IDS

				double parsed = Double.parseDouble(numbers[i]);
				nums.add(parsed);
			}

			MGFtable.add(nums);

		}
		sc.close();
		return MGFtable;

	}

	/**
	 * To construct an array of be.ugent.psb.setas.bdmodel.test branch lengths for the WGD event on the
	 * branch to the target child node
	 */
	public double[] setUpWGDtestBrLenxs(Node child) {

		double childBlen = child.getbLen();
		double[] WGDtestBrLenxs = new double[numPartitionsOnSusBr];

		/* brLen =0 corresponds to the parent node and brLen= brLen corresponds to the child node,
		 *  so we are not interested in the first and last elements of this array */
		
		double lenx_1Partition = childBlen / numPartitionsOnSusBr;
		
		for (int i = 1; i < WGDtestBrLenxs.length; i++) {

			WGDtestBrLenxs[i] = i * lenx_1Partition;
		}

		return WGDtestBrLenxs;
	}

	/**
	 * Calculates log-likelihood array for different Branch lengths of a WGD on
	 * ONE branch and for ONE MGF data
	 */
	public double[] calcLk_testBrLenxs_1MGF(int[] mgfData,
			double lambdaStar, int rootSizeStar) {

		/* first we build the array of WGD be.ugent.psb.setas.bdmodel.test-branchLengths */
		double[] WGDtestBrLenxs = setUpWGDtestBrLenxs(childToWGD);

		/*
		 * The corresponding Likelihoods to each of these positions are stored
		 * here:
		 */
		double[] lkForWGDWithTestBrLenx = new double[WGDtestBrLenxs.length];
		double[] logLkForWGDWithTestBrLenx = new double[WGDtestBrLenxs.length];

		String parentName = parentToWGD.getName();
		String childName = childToWGD.getName();

		Node wg = new Node();
		wg.isWGM = true;
		wg.multFactor = 2; /* For now we only be.ugent.psb.setas.bdmodel.test for WGDs */

		LikeLihood lk = new LikeLihood(lambdaStar,
				rootOfTree.getmaxNodeSize() + 1, probCache);

		// at index zero the wgd-blen is zero, and in the last index, it falls
		// on the child node, so not interested in those slots of the array. 
		for (int i = 1; i < WGDtestBrLenxs.length; i++) {

			wg.setbLen(WGDtestBrLenxs[i]);

			rootOfTree.addWGM(parentName, childName, wg);
			/* set the leaf values of the tree to that of the marker gene family */
			rootOfTree.setLeafValues(mgfData);

			double [] temp = lk.calcInternalLk(arln);
			 
			lkForWGDWithTestBrLenx[i] = temp[rootSizeStar];
			 
			logLkForWGDWithTestBrLenx[i] = Math.log(lkForWGDWithTestBrLenx[i]);

			rootOfTree = Node.resetAllSettings(treeFile, rootOfTree.getmaxNodeSize());
			rootOfTree.insertWGMsToTheTree(fixedWGDList);

		}
		return logLkForWGDWithTestBrLenx;

	}

	/**
	 * combines the resulting likelihoods arrays for different MGF data into one
	 * array by adding up log-liks
	 */

	public double[] calcCombLgLk_AllMGF_1Branch() {

		double[] logLk_1MGF_TestBrLenx = new double[numPartitionsOnSusBr];
		double[] sumLogLk_allMGF_TestBrLenx = new double[numPartitionsOnSusBr];

		/* read the MGF data file */
		List<List<Double>> MGF_data = readMGFdata();

		/* loop over MGF data */
		for (int mgf = 0; mgf < MGF_data.size(); mgf++) {

			rootOfTree = Node.resetAllSettings(treeFile, rootOfTree.getmaxNodeSize());
			rootOfTree.insertWGMsToTheTree(fixedWGDList);

			List<Double> mgfRecord = MGF_data.get(mgf);
			/* set the root size, lambda and data according to a specific MGF */
			int rootStar = mgfRecord.get(0).intValue();
			double lambdaStar = mgfRecord.get(1);
			
			int[] mgfGeneCounts = new int[mgfRecord.size() - 2];
			for (int m = 2; m < mgfRecord.size(); m++) {
				mgfGeneCounts[m - 2] = mgfRecord.get(m).intValue();
			}

			logLk_1MGF_TestBrLenx = calcLk_testBrLenxs_1MGF(
					mgfGeneCounts, lambdaStar, rootStar);

			for (int j = 1; j < numPartitionsOnSusBr; j++) {

				sumLogLk_allMGF_TestBrLenx[j] += logLk_1MGF_TestBrLenx[j];

			}
		}

		return sumLogLk_allMGF_TestBrLenx;
	}

	/* MLP Most Likely Place*/
	public double findSet_MLP_1WGD_1Branch_AllMGF(Branch b) {

		double[] lgLikelihoods = calcCombLgLk_AllMGF_1Branch();
		
		int maxIndex = FindMaxArray.findIndexOfMaxValue(lgLikelihoods);
		double bestWGDBranchLenght = maxIndex* (b.getBranchLenght() / numPartitionsOnSusBr);

		b.setMostLikelyPlaceWGD(bestWGDBranchLenght);
		return bestWGDBranchLenght;
	}

	/**********************/
	/** for 2 WGDs on One branch */
	public double[][] calcLk_testBrLenxs_1MGF_2WGDs(
			int[] markerFamilyData, double lambdaStar, int rootSizeStar) {

		/* first we build the array of WGD be.ugent.psb.setas.bdmodel.test-branchLengths */
		double[] WGDtestBrLenxs = setUpWGDtestBrLenxs(childToWGD);

		/*
		 * The corresponding Likelihoods to each of these positions are stored
		 * here:
		 */
		double[][] lkForWGDWithTestBrLenx = new double[WGDtestBrLenxs.length][WGDtestBrLenxs.length];
		double[][] logLkForWGDWithTestBrLenx = new double[WGDtestBrLenxs.length][WGDtestBrLenxs.length];

		String parentName = parentToWGD.getName();
		String childName = childToWGD.getName();

		Node wg1 = new Node();
		wg1.isWGM = true;
		wg1.multFactor = 2; /* For now we only be.ugent.psb.setas.bdmodel.test for WGDs */

		Node wg2 = new Node();
		wg2.isWGM = true;
		wg2.multFactor = 2;

		LikeLihood lk = new LikeLihood(lambdaStar,
				rootOfTree.getmaxNodeSize() + 1);

		/* at index zero the wgd-blen is zero, and in the last index, it falls
		 on the child node, so not interested in those slots of the array
		 It will be above-triagonal becasue of different combinations of 2
		 wgds on one branch */
		for (int i = 1; i < WGDtestBrLenxs.length; i++) {

			wg1.setbLen(WGDtestBrLenxs[i]);
			rootOfTree.addWGM(parentName, childName, wg1);
			
			/* To avoid placing two WGDs at the same branch length, 
			and avoid repeptetive symmetrical combinations */
			for (int j = i + 1; j < WGDtestBrLenxs.length; j++) {
				
				wg2.setbLen(WGDtestBrLenxs[j]-wg1.getbLen());
				rootOfTree.addWGM(wg1.getName(), childName, wg2);
				/*
				 * set the leaf values of the tree to that of the marker gene
				 * family
				 */
				rootOfTree.setLeafValues(markerFamilyData);

				 double [] lks= lk.calcInternalLk(arln);
				 
				lkForWGDWithTestBrLenx[i][j]= lks[rootSizeStar];
				logLkForWGDWithTestBrLenx[i][j] = Math
						.log10(lkForWGDWithTestBrLenx[i][j]);

				rootOfTree = Node.resetAllSettings(treeFile, rootOfTree.getmaxNodeSize());
				rootOfTree.insertWGMsToTheTree(fixedWGDList);

			}
		}
		return logLkForWGDWithTestBrLenx;

	}

	public double[][] calcComLk_AllMGF_testBrLenxs_2WGDs() {

		double[][] logLkForWGDsPositionedAtTestBrLenx = new double[numPartitionsOnSusBr][numPartitionsOnSusBr];
		double[][] sumLogLkForWGDsPositionedAtTestBrLenx = new double[numPartitionsOnSusBr][numPartitionsOnSusBr];

		/* read the MGF data file */
		List<List<Double>> MGF_data = readMGFdata();

		/* loop over MGF data */
		for (int mgfCount = 0; mgfCount < MGF_data.size(); mgfCount++) {

			rootOfTree = Node.resetAllSettings(treeFile, rootOfTree.getmaxNodeSize());
			rootOfTree.insertWGMsToTheTree(fixedWGDList);

			List<Double> mgfRecord = MGF_data.get(mgfCount);

			int[] mgfData = new int[mgfRecord.size() - 2];

			/* set the root size, lambda and data according to a specific MGF */
			int rootStar = mgfRecord.get(0).intValue();
			double lambdaStar = mgfRecord.get(1);

			for (int m = 2; m < mgfRecord.size(); m++) {
				mgfData[m - 2] = mgfRecord.get(m).intValue();
			}

			logLkForWGDsPositionedAtTestBrLenx = calcLk_testBrLenxs_1MGF_2WGDs(
					mgfData, lambdaStar, rootStar);

			for (int i = 1; i < numPartitionsOnSusBr; i++) {
				for (int j = i+1; j < numPartitionsOnSusBr; j++) {

					// add up log-lks from diffrent mgf data for one combination
					// of (i,j) wgd on the branch
					sumLogLkForWGDsPositionedAtTestBrLenx[i][j] += logLkForWGDsPositionedAtTestBrLenx[i][j];
				}

			}
		}

		return sumLogLkForWGDsPositionedAtTestBrLenx;
	}

	/* returns the most likely branch lengths of 2 WGDs on once branch accroding to all MGF data,
	 * IN THE WAY IT SHOULD BE IN WGD FILE: 2nd WGD BrLen = distance to the first(older) one */
	public double[] findSet_MLPs_2WGD_1Branch(Branch b) {

		double[][] likelihoods = calcComLk_AllMGF_testBrLenxs_2WGDs();

		int[] maxIndexes = new int[2];
		/* returns (i,j) corresponding to maxLogLk(i,j)*/
		maxIndexes = FindMaxArray.findIndexOfMaxValueInMatrix(likelihoods,
				numPartitionsOnSusBr, numPartitionsOnSusBr);

		double[] bestWGDBrLenx = new double[2];
		bestWGDBrLenx[0] = maxIndexes[0]
				* (b.getBranchLenght() / numPartitionsOnSusBr);

		bestWGDBrLenx[1] = maxIndexes[1]
				* (b.getBranchLenght() / numPartitionsOnSusBr)- bestWGDBrLenx[0];

		b.setMostLikelyPlacesOf2WGDs(bestWGDBrLenx);

		return bestWGDBrLenx;
	}

	// method: determine if a branch should have one or two WGDs

	/**********************/

	/**
	 * calculates the combined likelihood of all MGF data, for a fixed tree with
	 * WGDs on it: To use in the next method!
	 */

	public double calcComLkAllMGF_tree_AddedNewWGDs() {

		List<List<Double>> MGF_data = readMGFdata();
		double sumLogLk = 0;

		for (int i = 0; i < MGF_data.size(); i++) {

			List<Double> mgfRecord = MGF_data.get(i);
			int[] mgfGeneCounts = new int[mgfRecord.size() - 2];

			int rootStar = mgfRecord.get(0).intValue();
			// System.out.println("rootStar " + rootStar);

			double lambdaStar = mgfRecord.get(1);
			// System.out.println("lambdaStar" + lambdaStar);

			for (int m = 2; m < mgfRecord.size(); m++) {
				mgfGeneCounts[m - 2] = mgfRecord.get(m).intValue();
			}

			LikeLihood lk = new LikeLihood(lambdaStar,
					rootOfTree.getmaxNodeSize());

			rootOfTree.setValue(rootStar);
			rootOfTree.setLeafValues(mgfGeneCounts);

			double [] lks = lk.calcInternalLk(arln);
			sumLogLk += lks[rootStar];
					
		}

		return sumLogLk;
	}

	/*
	 * Calculate different likelihoods for one MGF data by positioning WGD on
	 * different sub-branches of the tree, always placing the WGD in the middle
	 * of the branch or calculate the most likely position first
	 */

	public ArrayList<BranchList_Likelihood> calcMaxLk_KOutNbranchAllMGF(int k) {

		ArrayList<BranchList_Likelihood> brList_Lk = new ArrayList<BranchList_Likelihood>();

		ArrayList<Branch> allBranches = rootOfTree.findAllBranches(rootOfTree);

		ChooseKoutOfN<Branch> choice_KoutN_br = new ChooseKoutOfN<Branch>(allBranches);
		List<List<Branch>> subsets_kBrs = choice_KoutN_br.generateSubSets(k);

		// for(List<Branch> brl : subsetOfKBranches){
		// System.out.println("branch list: ");
		// for(Branch brc : brl){
		// System.out.println("parent: "+ brc.getParent().getName()+ " child: "+
		// brc.getChild().getName());}
		// }

		int sizeOfKCombinations = subsets_kBrs.size();

		double[] lkOfEachCombination = new double[sizeOfKCombinations];

		for (int j = 0; j < sizeOfKCombinations; j++) {

			List<Branch> kBranches = subsets_kBrs.get(j);

			// add all the WGDs on their most likely place on the branches
			for (int br = 0; br < k; br++) {

				Branch br_k = kBranches.get(br);
				
				double mlp_k = findSet_MLP_1WGD_1Branch_AllMGF(br_k);
				
				// br_k.setMostLikelyPlaceWGD(mlp_k);
				rootOfTree.addWGMOnaSpecificBranch(br_k, mlp_k);
			}

			lkOfEachCombination[j] = calcComLkAllMGF_tree_AddedNewWGDs();
			
			//Here we have to reset all the changes to the tree before we move to the next combination of k branches

		}
		int maxIndex = FindMaxArray.findIndexOfMaxValue(lkOfEachCombination);

		List<Branch> bestBrList = subsets_kBrs.get(maxIndex);
		double correspondingLk = lkOfEachCombination[maxIndex];

		brList_Lk.add(new BranchList_Likelihood(bestBrList, correspondingLk));

		return brList_Lk;

	}

	private class BranchList_Likelihood {
		public List<Branch> branchList;
		public double lk;

		public BranchList_Likelihood(List<Branch> branchList, double lk) {
			this.branchList = branchList;
			this.lk = lk;
		}
	}

	/*
	 * Finds optimal branch-likelihood list for any combination of k <
	 * numbBranchInTheTree and returns the most likely branch-lk list
	 */
	public ArrayList<BranchList_Likelihood> calcMaxLkCombWGDs() {

		ArrayList<Branch> arb = rootOfTree.findAllBranches(rootOfTree);
		int numBranches = arb.size();

		double maxLk = 0;
		int maxIndex = 0;

		for (int k = 1; k < numBranches; k++) {

			ArrayList<BranchList_Likelihood> br_lk = calcMaxLk_KOutNbranchAllMGF(k);

			if (br_lk.get(0).lk > maxLk) {
				maxLk = br_lk.get(0).lk;
				maxIndex = k;
			}

		}

		return calcMaxLk_KOutNbranchAllMGF(maxIndex);
	}

	public static void main(String[] args) {

		long startTime = System.currentTimeMillis();

		double partitionSize = 0.1;
		int defaultmaxNodeSize = 100;
		
		Node root = SpeciesTreeParser.buildInsertWGDsandPartitionTree(args[0], args[1],
				partitionSize, defaultmaxNodeSize);

		ReadGFcountsFile rgf = new ReadGFcountsFile();
		List<List<Integer>> nonRepetetiveCounts = rgf
				.read_unique(args[2]);
		
		ArrayList<String> gfIDs = rgf.getGfIDs_unique();

		ArrayList<Node> leaves = root.getLeaves();

		// for(int i=0; i<leaves.size();i++){
		// System.out.println(leaves.get(i).getName());
		// }

		// System.out.println(leaves.size());
		
		/* If MGF1 has maxCount = 100 and MGF2 has it = 200, when combining the lk for slot 101 in lkArray at an internal node, MGF1 will multiply lk by 0, 
		 * resulting in eliminating this slot form the rest of the competition. So we should set the maxNodeSize to Min between all MGFs.
		 * In reality, such hige GFs will not be MGF, so not a big concern.*/
		int mgf =0; 
	    SpeciesTreeParser.setLeavesValues(root, nonRepetetiveCounts, mgf);
		ArrayList<Node> speciesTree = SpeciesTreeParser.setMaxNodeSizeAccToGF(root, nonRepetetiveCounts, mgf);
		
		
		ProbCalculator probCache = new ProbCalculator();

		ReverseEngWGDs rvsWGD = new ReverseEngWGDs(args[0], args[2], root, probCache);

		rvsWGD.numPartitionsOnSusBr = 2;
		rvsWGD.parentToWGD = root.findNodeWithName("4");
		// System.out.println("parent blen: "+rvsWGD.parentToWGD.getbLen());

		rvsWGD.childToWGD = root.findNodeWithName("Brap");
		// System.out.println("childblen: " + rvsWGD.childToWGD.getbLen());

		// double[] log_lks = rvsWGD.calcComLkAllMGFAllPosWGDOn1Branch();

		// for (int i = 0; i < log_lks.length; i++) {
		// System.out.println("log_lks: "+log_lks[i]);
		// }

		ArrayList<BranchList_Likelihood> arl_brLk = rvsWGD
				.calcMaxLk_KOutNbranchAllMGF(8);

		for (int i = 0; i < arl_brLk.size(); i++) {

			for (int j = 0; j < arl_brLk.get(i).branchList.size(); j++) {

				System.out.println("lk: "
						+ arl_brLk.get(i).lk
						+ " parent: "
						+ arl_brLk.get(i).branchList.get(j).getParent()
								.getName()
						+ arl_brLk.get(i).branchList.get(j).getChild()
								.getName());
			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("\t" + (endTime - startTime));
	}

}
