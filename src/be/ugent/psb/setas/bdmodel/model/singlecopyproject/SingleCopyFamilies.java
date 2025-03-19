package be.ugent.psb.setas.bdmodel.model.singlecopyproject;
//
//import java.io.FileNotFoundException;
//import java.io.FileReader;
//import java.util.ArrayList;
//import java.util.List;
//import java.util.Map;
//import java.util.Random;
//import java.util.Scanner;
//import java.util.TreeMap;
//
//import be.ugent.psb.setas.bdmodel.model.GenerateObservations;
//import be.ugent.psb.setas.bdmodel.model.Node;
//import be.ugent.psb.setas.bdmodel.model.ProbCalculator;
//import be.ugent.psb.setas.bdmodel.parsers.NewickParser;
//import be.ugent.psb.setas.bdmodel.parsers.WGMparser;
//
public class SingleCopyFamilies {
//
//	public SingleCopyFamilies(Node root, int numObsPerLam, int testSize,
//			int numOfBins, int maxNumberOfZeros) {
//		super();
//		this.root = root;
//		this.numObsPerLam = numObsPerLam;
//		this.testSize = testSize;
//		this.numOfBins = numOfBins;
//		this.maxNumberOfZeros = maxNumberOfZeros;
//	}
//
//	private Node root;
//
//	private int numObsPerLam;
//	private int testSize;
//	private int numOfBins;
//	private int maxNumberOfZeros;
//
//	public Node getRoot() {
//		return root;
//	}
//
//	public void setRoot(Node root) {
//		this.root = root;
//	}
//
//	public int getNumObsPerLam() {
//		return numObsPerLam;
//	}
//
//	public void setNumObsPerLam(int numObsPerLam) {
//		this.numObsPerLam = numObsPerLam;
//	}
//
//	public int getTestSize() {
//		return testSize;
//	}
//
//	public void setTestSize(int testSize) {
//		this.testSize = testSize;
//	}
//
//	public int getNumOfBins() {
//		return numOfBins;
//	}
//
//	public void setNumOfBins(int numOfBins) {
//		this.numOfBins = numOfBins;
//	}
//
//	public int getMaxNumberOfZeros() {
//		return maxNumberOfZeros;
//	}
//
//	public void setMaxNumberOfZeros(int maxNumberOfZeros) {
//		this.maxNumberOfZeros = maxNumberOfZeros;
//	}
//
//	public double[] readSampleLambdas(String fileName) {
//		FileReader fin = null;
//		try {
//			fin = new FileReader(fileName);
//		} catch (FileNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//
//		Scanner sc = new Scanner(fin);
//		sc.nextLine(); /* The first line is a header: */
//
//		List<Double> lambdas = new ArrayList<Double>();
//
//		while (sc.hasNextLine()) {
//			String line = sc.nextLine();
//
//			if (!line.equalsIgnoreCase("")) {
//
//				double parsed = Double.parseDouble(line);
//
//				lambdas.add(parsed);
//			}
//		}
//		sc.close();
//
//		double[] lambdasArray = new double[lambdas.size()];
//
//		for (int i = 0; i < lambdas.size(); i++) {
//
//			lambdasArray[i] = lambdas.get(i);
//		}
//		return lambdasArray;
//	}
//
//	public List<List<Double>> readTabDelimitedFile(String filename) {
//
//		FileReader fin = null;
//		try {
//			fin = new FileReader(filename);
//		} catch (FileNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//
//		Scanner sc = new Scanner(fin);
//
//		List<List<Double>> table = new ArrayList<List<Double>>();
//
//		while (sc.hasNextLine()) {
//
//			String line = sc.nextLine();
//			String[] numbers = line.split("\t");
//
//			/* Now we have the numbers as strings in the array "numbers", so: */
//			List<Double> nums = new ArrayList<Double>();
//
//			for (int i = 0; i < numbers.length; i++) {
//
//				double parsed = Double.parseDouble(numbers[i]);
//				nums.add(parsed);
//			}
//			table.add(nums);
//
//		}
//		sc.close();
//		return table;
//	}
//
//	public List<List<Double>> reversetable(List<List<Double>> tab) {
//
//		int length = tab.get(0).size();
//		List<List<Double>> tableReverse = new ArrayList<List<Double>>();
//
//		for (int j = 0; j < length; j++) {
//			tableReverse.add(new ArrayList<Double>());
//		}
//
//		for (int i = 0; i < tab.size(); i++) {
//			for (int j = 0; j < tab.get(i).size(); j++) {
//				List<Double> row = tableReverse.get(j);
//				row.add(tab.get(i).get(j));
//			}
//		}
//		return tableReverse;
//	}
//
//	/* Takes an observation array and return SC percentage */
//	public double calcSCpercent(int[] observations) {
//
//		double percentage = 0;
//		int counter = 0;
//
//		int numberOfZeros = 0;
//
//		for (int i = 0; i < observations.length; i++) {
//
//			if (observations[i] == 1) {
//
//				counter += 1;
//			}
//
//			if (observations[i] == 0) {
//
//				numberOfZeros += 1;
//			}
//		}
//
//		int realObservationLenx = observations.length - numberOfZeros;
//
//		/* To have percentages as values between 0 and 100 */
//
//		percentage = (counter * numOfBins) / (realObservationLenx);
//		System.out.println(percentage);
//
//		return percentage;
//
//	}
//
//	/*
//	 * Takes an array of calculated SC percentages and creats percentages bins
//	 * and fills in them
//	 */
//	public int[] creatBins(double[] singleCopyPercentage) {
//
//		int[] bins = new int[numOfBins];
//
//		for (int j = 0; j < singleCopyPercentage.length; j++) {
//
//			/*
//			 * for the analysis of zhen the sc-percentages are between 0 and 1
//			 * so i multiply by 100, for my analysis shoudl be removed
//			 */
//			int indexInBins = (int) Math.floor(singleCopyPercentage[j]);
//
//			/* if you reached the end value (100), add it to previous bin (99): */
//			if (indexInBins == numOfBins) {
//				bins[indexInBins - 1] += 1;
//			} else {
//				bins[indexInBins] += 1;
//			}
//		}
//
//		return bins;
//	}
//
//	/*
//	 * combines the above two methods: generates observation for one lambda,
//	 * saves SC percentages, passes it to second method above, to calculate bins
//	 */
//	public int[] creatBinForOneLam(double lambda, ProbCalculator probCache) {
//
//		int[] bins = new int[numOfBins];
//
//		GenerateObservations go = new GenerateObservations(0,
//				root.getmaxNodeSize(), false, probCache);
//
//		double[] scPercentages = new double[numObsPerLam];
//
//		for (int counter = 0; counter < numObsPerLam; counter++) {
//
//			double[] qWGM = genRandRetentionRate();
//
//			// System.out.println("qWGD: "+qWGM[0]+"  qWGT: "+qWGM[1]);
//
//			int[] obsv = go.generateObservation_SC(root, testSize, lambda,
//					maxNumberOfZeros, qWGM[0], qWGM[1]);
//
//			scPercentages[counter] = calcSCpercent(obsv);
//
//		}
//
//		bins = creatBins(scPercentages);
//		return bins;
//	}
//
//	public int[] sumOfTwoBins(int[] binOne, int[] binTwo) {
//
//		if (binOne.length != binTwo.length) {
//			System.out
//					.println("Warning: adding up two sets of bins with unequal lenght. Only the length of the fist bin set is used.");
//		}
//		int[] binSum = new int[binOne.length];
//
//		for (int i = 0; i < binOne.length; i++) {
//
//			binSum[i] = binOne[i] + binTwo[i];
//
//		}
//		return binSum;
//	}
//
//	public int[] creatFinalBins(double[] lambdas, ProbCalculator probCache) {
//		int[] sumBins = creatBinForOneLam(lambdas[0], probCache);
//		
//		for (int lam = 1; lam < lambdas.length; lam++) {
//
//			int[] newBin = creatBinForOneLam(lambdas[lam], probCache);
//
//			sumBins = sumOfTwoBins(sumBins, newBin);
//
//		}
//
//		return sumBins;
//	}
//
//	public double[] creatUniformSampleLambdas(double minLambda,
//			double maxLambda, int numberOfSamplePoints) {
//
//		double[] lambdas = new double[numberOfSamplePoints];
//		Random random = new Random();
//
//		for (int i = 0; i < lambdas.length; i++) {
//			/*
//			 * double randomValue = rangeMin + (rangeMax - rangeMin) *
//			 * r.nextDouble()
//			 */
//			lambdas[i] = minLambda + (maxLambda - minLambda)
//					* random.nextDouble();
//		}
//
//		return lambdas;
//
//	}
//
//	public double[] genRandRetentionRate() {
//
//		Random random = new Random();
//
//		double[] ratesWGM = new double[2];
//
//		ratesWGM[0] = 1 + random.nextDouble();
//		ratesWGM[1] = 1.5 * ratesWGM[0];
//
//		return ratesWGM;
//	}
//
////	public static void main(String args[]) {
////
////		int maxNodeSize = 100;
////
////		NewickParser np = new NewickParser();
////		Node root = np.buildAndPartitionTree(args[0], maxNodeSize);
////
////		WGMparser wgm = new WGMparser();
////		List<List<String>> wgdList = wgm.readInputFile(args[1]);
////		root.insertWGMsToTheTree(wgdList);
////
////		ArrayList<Node> leaves = root.getLeaves();
////		root.setLeaves(leaves);
////		root.setNumberOfLeaves(leaves.size());
////
////		// Queue<Double> q = root.allBlen();
////		// int numGeneFamilies = 9178;
////
////		SingleCopyFamilies scf = new SingleCopyFamilies(root, 1, 1, 100, 5);
////		// double [] lambdas = scf.creatUniformSampleLambdas(0, 9.99,
////		// numGeneFamilies);
////
////		// double [] lambdas = scf.readSampleLambdas(args[2]);
////		List<List<Double>> table = scf.readTabDelimitedFile(args[2]);
////		List<List<Double>> lambdasTable = scf.reversetable(table);
////		List<List<Double>> scpDistributions = scf.reversetable(table);
////
////		// int[] sumBins = new int [100];
////
////		for (int i = 0; i < 1000; i++) {
////
////			// System.out.println(i);
////			List<Double> li = scpDistributions.get(i);
////			double[] scps = new double[li.size()];
////
////			for (int m = 0; m < li.size(); m++) {
////				scps[m] = li.get(m);
////			}
////
////			int[] scpsNewBin = scf.creatBins(scps);
////
////			for (int bin : scpsNewBin) {
////				System.out.print(bin + "\t");
////				// System.out.print("\t");
////			}
////			System.out.println();
////
////			// sumBins = scf.sumOfTwoBins(sumBins, scpsNewBin);
////
////		}
////
////		// for(int i=0; i<sumBins.length;i++){
////		// System.out.println(sumBins[i]/1000);}
////
////		int sampleNom = Integer.parseInt(args[3]);
////		List<Double> li = lambdasTable.get(sampleNom);
////		double[] lambdas = new double[li.size()];
////		for (int m = 0; m < li.size(); m++) {
////			lambdas[m] = li.get(m);
////
////		}
////
////		// GenerateBiasedSampleOfLambdas gbsl = new
////		// GenerateBiasedSampleOfLambdas();
////		// gbsl.sampleSize = 9178; GENERATE THIS MANY FAMILIES: DO THE CALCS ON
////		// STEPSIZE ETC
////		// double[] lambdas2 = gbsl.generateSampleLambdas(args[1]);
////		Map<CacheKey, Double> probCache = new TreeMap<ProbCalculator.CacheKey, Double>();
////		int[] sumBins = scf.creatFinalBins(lambdas, probCache);
////
////		// double [] Zhen_sc_percentages =
////		// scf.readSampleLambdas("/home/setas/workspace/BirthDeathModel/src/Files/Zhen_scPercentages.txt");
////		// int [] bins = scf.creatBins(Zhen_sc_percentages);
////
////	}
//
}
