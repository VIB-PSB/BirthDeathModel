package be.ugent.psb.setas.bdmodel.model.RVS_Engineering;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ugent.psb.setas.bdmodel.model.Branch;
import be.ugent.psb.setas.bdmodel.model.LikeLihood;
import be.ugent.psb.setas.bdmodel.model.MathOperations;
import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.bdmodel.model.ProbCalculator;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;

public class ReverseEngWGDs_Simple {

	private Node rootOfTree;

	private String treeFile;
	private String WGDsfile;
	private String mgfFile;

	private int ignoreLineNumWGDFile;

	private ProbCalculator probCache;
	private int maxNodeSize;

	// To include all lines of WGD file:
	public ReverseEngWGDs_Simple(String treeFile, String WGDsfile,
			String mgfFile, ProbCalculator probCache, int maxNodeSize,
			double partitionLenx) {

		this.treeFile = treeFile;
		this.WGDsfile = WGDsfile;
		this.mgfFile = mgfFile;
		this.maxNodeSize = maxNodeSize;
		this.rootOfTree = SpeciesTreeParser.buildInsertWGDsandPartitionTree(treeFile,
				WGDsfile, partitionLenx, maxNodeSize);
		this.probCache = probCache;
	}

	// To remove one line from the WGD file:
	public ReverseEngWGDs_Simple(String treeFile, String WGDsfile,
			String mgfFile, int ignoreLine, ProbCalculator probCache,
			int maxNodeSize, double partitionLenx) {

		this.treeFile = treeFile;
		this.WGDsfile = WGDsfile;
		this.mgfFile = mgfFile;
		this.maxNodeSize = maxNodeSize;
		this.ignoreLineNumWGDFile = ignoreLine;
		this.rootOfTree = SpeciesTreeParser.buildAndPartitionTree_ReverseEng(
				treeFile, WGDsfile, partitionLenx, maxNodeSize,
				ignoreLineNumWGDFile);
		this.probCache = probCache;
	}

	/**
	 * The MGF file should be a file in each line one Marker-Gene-Family ID and
	 * tab separated: rootSize*, lambda*, and then the gene counts
	 */
	public List<List<Double>> readMGFgeneCounts() {

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

			for (int i = 1; i < numbers.length; i++) { // i=0 GF-IDS i=1
														// rootSize i=2 lambda
														// i=3 loglks i=4 pv
														// i=5,6.. gene counts

				double parsed = Double.parseDouble(numbers[i]);
				nums.add(parsed);
			}

			MGFtable.add(nums);

		}
		sc.close();
		return MGFtable;

	}

	public double calcLogLk_1MGF(List<Double> mgfRecord, ArrayList<Node> arln) {

		int rootStar = mgfRecord.get(0).intValue();
		double lambdaStar = mgfRecord.get(1);

		int[] mgfGeneCounts = new int[mgfRecord.size() - 4];
		for (int m = 4; m < mgfRecord.size(); m++) { // 0=r* 1=lam* 2=lglk*
														// 3=pv*
			mgfGeneCounts[m - 4] = mgfRecord.get(m).intValue();

		}

		rootOfTree.setLeafValues(mgfGeneCounts);

		LikeLihood lkhd = new LikeLihood(lambdaStar,
				rootOfTree.getmaxNodeSize() + 1, probCache);

		double[] Lk_1MGF = lkhd.calcInternalLk(arln);

		double[] logLk_1MGF = MathOperations.giveLogArray(Lk_1MGF);

		return logLk_1MGF[rootStar];

	}

	public double calcCombLgLk_AllMGF(Node root) {

		double sumLogLk_allMGF = 0;

		List<List<Double>> MGF_data = readMGFgeneCounts();
		ArrayList<Node> arln = SpeciesTreeParser.setMaxNodeSize(root,
				this.maxNodeSize);

		for (int mgf = 0; mgf < MGF_data.size(); mgf++) {

			List<Double> mgfRecord = MGF_data.get(mgf);
			double logLk_1MGF = calcLogLk_1MGF(mgfRecord, arln);

			sumLogLk_allMGF = sumLogLk_allMGF + logLk_1MGF;
		}

		return sumLogLk_allMGF;
	}

	/*******************************/
	// To add a WGD on one branch:

	public static int findIndexOfMaxValue_editted(double[] a) {
		int index = 0;
		for (int i = 1; i < a.length; i++) {
			if (a[i] > a[index]) {
				index = i;
			}
		}
		return index;
	}

	public void calcMaxLk_addWGD_1branchAllMGF() {

		ArrayList<Branch> allBranches = rootOfTree.findAllBranches(rootOfTree);

		double[] save_lks = new double[allBranches.size()];

		for (int j = 0; j < allBranches.size(); j++) {

			Branch br = allBranches.get(j);

			// add all the WGDs on the middle of a branches
			rootOfTree.addWGMOnaSpecificBranch(br, br.getBranchLenght() / 2);

			for (Node n : rootOfTree.postOrder()) {
				if (n.isWGD) {
					System.out.println(n.getName());
				}
			}
			// I changed this function to be.ugent.psb.setas.bdmodel.test removing WGDs for MGFs one by one
			save_lks[j] = calcCombLgLk_AllMGF(rootOfTree);

			System.out.println(save_lks[j]);

			// Here we have to reset all the changes to the tree before we move
			// to the next combination of k branches
			rootOfTree = Node.resetAllSettings_new(treeFile, WGDsfile, 0.1,
					this.maxNodeSize);
		}
		int maxIndex = findIndexOfMaxValue_editted(save_lks);

		Branch bestBr = allBranches.get(maxIndex);
		double correspondingLk = save_lks[maxIndex];

//		System.out.println(correspondingLk);
//
//		System.out.println(bestBr.getParent().getName() + "  "
//				+ bestBr.getChild().getName());

	}

	public static void main(String[] args) {

		ProbCalculator probCache = new ProbCalculator();
		int maxNodeSize = 100;
		double partitioningSize = 0.1;

		// ReverseEngWGDs_Simple rvsEng_original = new
		// ReverseEngWGDs_Simple(args[0],args[1],args[2],probCache);
		//
		// double sumLkAllMGF =
		// rvsEng_original.calcCombLgLk_AllMGF(rvsEng_original.rootOfTree);
		//
		// System.out.println(sumLkAllMGF);

		// for(int mgf=0; mgf <709; mgf++){

		for (int i = 0; i <= 17; i++) { // because we have 17 lines in our WGD
										// file = 18 WGD/T s

			if (i == 8 || i == 11 || i == 15) { // older WGD require another WGD
												// file to be used: args[2]
												// instead of args[1]

				ReverseEngWGDs_Simple rvsEng2 = new ReverseEngWGDs_Simple(
						args[0], args[2], args[3], i, probCache, maxNodeSize,
						partitioningSize);

				double sumLkAllMGF_testWGDset2 = rvsEng2
						.calcCombLgLk_AllMGF(rvsEng2.rootOfTree);
				System.out.println(sumLkAllMGF_testWGDset2);
			}
			
			else{
			// Removing WGD at line i in WGDfile:
			ReverseEngWGDs_Simple rvsEng = new ReverseEngWGDs_Simple(args[0],
					args[1], args[3], i, probCache, maxNodeSize,
					partitioningSize);
			// Test 1: calculate combined lk
			double sumLkAllMGF_testWGDset = rvsEng
					.calcCombLgLk_AllMGF(rvsEng.rootOfTree);
			System.out.println(sumLkAllMGF_testWGDset);
			}

		}

		// Test 2: try to put one WGD on a branch, take the branch with max Lk:
		// rvsEng.calcMaxLk_addWGD_1branchAllMGF();

		// }

	}
}
