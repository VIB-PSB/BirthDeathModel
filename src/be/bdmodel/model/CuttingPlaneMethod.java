package be.ugent.psb.setas.bdmodel.model;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Queue;
import java.util.Random;
import java.util.Scanner;

import be.ugent.psb.setas.bdmodel.parsers.ReadGFcountsFile;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;
import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class CuttingPlaneMethod {

	private List<Node> aryln;
	private LikeLihood lk;
	private int rootSize;
	private List<List<Integer>> gfCounts;
	private Node root;

	private final double stepSize;
	private final double deltaLocalMoves;
	private final double toleranceDerivative;
	private final double toleranceF;
	private final double precisionLambda;
	private final double minInterval;
	private final double maxInterval;
	private final double minAllowed;
	private final double maxAllowed;
	public int counter;
	private double optimalLambda;
	private double fOptimalLambda;
	private ProbCalculator probCalc;
	private HashMap<Double, Double> fCache;
	// private HashMap<Double, Double> fcombinedLoglkCache;

	public double getDeltaLocalMove() {
		return deltaLocalMoves;
	}

	public boolean isEqual(double a, double b, double PRECISION) {
		return (Math.abs(a - b) < PRECISION);
	}

	public int getRootSize() {
		return this.rootSize;
	}

	public double getOptimalLambda() {
		return optimalLambda;
	}

	public double getFoptimalLambda() {
		return fOptimalLambda;
	}

	public double getMaxInterval() {
		return this.maxInterval;
	}

	public double getMinInterval() {
		return minInterval;
	}

	public CuttingPlaneMethod(List<Node> speciesTree, int rootSize, double stepSize, double deltaLocalMoves,
			double tolD, double tolF, double minInterval, double maxInterval, double minAllowed, double maxAllowed,
			double precision, ProbCalculator probCache) {

		this.stepSize = stepSize;
		this.deltaLocalMoves = deltaLocalMoves;
		this.rootSize = rootSize;
		this.toleranceDerivative = tolD;
		this.toleranceF = tolF;
		this.minInterval = minInterval;
		this.maxInterval = maxInterval;
		this.minAllowed = minAllowed;
		this.maxAllowed = maxAllowed;
		this.precisionLambda = precision;
		this.probCalc = probCache;
		this.fCache = new HashMap();
		this.aryln = speciesTree;
	}

	public CuttingPlaneMethod(List<Node> speciesTree, int rootSize, double stepSize, double deltaLocalMoves,
			double tolD, double tolF, double minInterval, double maxInterval, double minAllowed, double maxAllowed,
			double precision, ProbCalculator probCache, Node root, List<List<Integer>> gfCounts) {

		this.stepSize = stepSize;
		this.deltaLocalMoves = deltaLocalMoves;
		this.rootSize = rootSize;
		this.toleranceDerivative = tolD;
		this.toleranceF = tolF;
		this.minInterval = minInterval;
		this.maxInterval = maxInterval;
		this.minAllowed = minAllowed;
		this.maxAllowed = maxAllowed;
		this.precisionLambda = precision;
		this.probCalc = probCache;
		this.fCache = new HashMap();
		// this.fcombinedLoglkCache= new HashMap();
		this.aryln = speciesTree;
		// this.root = root;
		// this.gfCounts = gfCounts;
	}

	public CuttingPlaneMethod(List<Node> speciesTree, int rootSize, double stepSize, double deltaLocalMoves,
			double tolD, double tolF, double minInterval, double maxInterval, double minAllowed, double maxAllowed,
			double precision, ProbCalculator probCache, Node root) {

		this.stepSize = stepSize;
		this.deltaLocalMoves = deltaLocalMoves;
		this.rootSize = rootSize;
		this.toleranceDerivative = tolD;
		this.toleranceF = tolF;
		this.minInterval = minInterval;
		this.maxInterval = maxInterval;
		this.minAllowed = minAllowed;
		this.maxAllowed = maxAllowed;
		this.precisionLambda = precision;
		this.probCalc = probCache;
		this.fCache = new HashMap();
		this.aryln = speciesTree;
		this.root = root;
	}

	// Original f
	public double f(double lambda) {

		Double value = fCache.get(lambda);
		if (value != null) {
			return value;
		}
		LikeLihood lk = new LikeLihood(lambda, aryln.get(0).getmaxNodeSize() + 1, probCalc);
		double[] lks = lk.calcInternalLk(aryln);

		double[] logLks = MathOperations.giveLogArray(lks);
		value = logLks[rootSize];
		fCache.put(lambda, value);
		return logLks[rootSize];
	}
	// To calculate the combined loglk
	// public double f(double lambda ) {
	// Double value = fcombinedLoglkCache.get(lambda);
	// if (value != null) {
	// return value;
	// }
	// LikeLihood lk = new LikeLihood(lambda,
	// aryln.get(0).getmaxNodeSize() + 1, probCalc);
	// double[] sumLogLks = lk.calcInternalLk_combLogLk(root, gfCounts, 0, 708,
	// aryln);
	//
	// value = sumLogLks[rootSize];
	//// System.out.println("lambda "+lambda+" sumLoglks: "+value);
	//
	// fCache.put(lambda, value);
	// return sumLogLks[rootSize];
	// }

	public double forwardDiff(double lambda) {
		double forwardDiff = 0;

		double fLambda = f(lambda);
		double fLambdaPlusDelta = f(lambda + stepSize);

		forwardDiff = (fLambdaPlusDelta - fLambda) / stepSize;

		return forwardDiff;
	}

	public double backwardDiff(double lambda) {
		double backwardDiff = 0;

		double fLambda = f(lambda);
		double fLambdaMinusDelta = f(lambda - stepSize);

		backwardDiff = (fLambda - fLambdaMinusDelta) / stepSize;

		return backwardDiff;
	}

	public double[] lineChars(double lambda, int fb) {

		/**
		 * lineChars[0] = lambda lineChars[1]= f(lambda) lineChars[2]= f'+ (lambda) or
		 * f'-(lambda)
		 */
		double[] lineChars = new double[3];

		lineChars[0] = lambda;
		lineChars[1] = f(lambda);
		// fb= 2 means calculate forward derivative and fb =3 means calculate backward
		// derivative
		if (fb == 2) {
			lineChars[2] = forwardDiff(lambda);
		} else if (fb == 3) {
			lineChars[2] = backwardDiff(lambda);
		}

		return lineChars;
	}

	public double findInterSectionLines(double[] line1, double[] line2) {
//		System.out.println("begining find Intersection  "+ line1[0]+"\t"+line2[0]);

		double lambdaStar = ((line1[2] * line1[0]) - (line2[2] * line2[0]) + line2[1] - line1[1])
				/ (line1[2] - line2[2]);

		double[] noNaN = solveNaN(lambdaStar, line1[0], line2[0]);
		return noNaN[0]; // lambdaStar_editted

//		System.out.println("end find Intersection  "+ lambdaStar+"\t"+line1[0]+"\t"+line2[0]);
//		return lambdaStar;
	}

	/* returns [lam*_1 , lam*_2] */
	public double[] calcLambdaStar12(double[] lineStar, double[] line1, double[] line2) {
//		System.out.println("begining of calc lambda star12 "+lineStar[0]+"\t"+line1[0]+"\t"+line2[0]);	
		double[] noNaN = solveNaN(lineStar[0], line1[0], line2[0]);

		lineStar = lineChars(noNaN[0], 2);
		line1 = lineChars(noNaN[1], 2);
		line2 = lineChars(noNaN[2], 2);

		double[] lambdaStar12 = new double[2];
		double lambdaStar_1 = findInterSectionLines(lineStar, line1);
		double lambdaStar_2 = findInterSectionLines(lineStar, line2);

		lambdaStar12[0] = lambdaStar_1;
		lambdaStar12[1] = lambdaStar_2;

//		System.out.println("end of calc lambda star12 "+lambdaStar12[0]+"\t"+lambdaStar12[1]);
		return lambdaStar12;
	}

	public double argMax(double lam1, double lam2, double lamStar) {

		if (Double.isNaN(lam1) || Double.isNaN(lam2) || Double.isNaN(lamStar)) {
			throw new RuntimeException("lambda = NaN");
		}

		final double fLam1 = f(lam1);
		final double fLamStar = f(lamStar);
		final double fLam2 = f(lam2);

		if (fLam1 <= fLamStar && fLam2 <= fLamStar) {
			return lamStar; // sub-optimality
		} else if (fLam1 <= fLamStar && fLam2 > fLamStar) {
			return lam2;
		}

		else if (fLam1 > fLamStar && fLam2 <= fLamStar) {
			return lam1;
		}

		else if (fLam1 > fLamStar && fLam2 > fLamStar) {
			if (fLam1 >= fLam2) {
				return lam1;
			} else {
				return lam2;
			}
		}
		System.err.println("argMax = -1 " + lam1 + "   " + lam2 + "   " + lamStar);
		return -1;
	}

	public int argMaxIndex(double lam1, double lam2, double lamStar) {

		if (Double.isNaN(lam1) || Double.isNaN(lam2) || Double.isNaN(lamStar)) {
			throw new RuntimeException("NaN is not good!");
		}

		final double fLam1 = f(lam1);
		final double fLamStar = f(lamStar);
		final double fLam2 = f(lam2);

		if (fLam1 <= fLamStar && fLam2 <= fLamStar) {
			return 2; // sub-optimality
		}

		else if (fLam1 <= fLamStar && fLam2 > fLamStar) {
			return 1;
		} else if (fLam1 > fLamStar && fLam2 <= fLamStar) {
			return 0;
		} else if (fLam1 > fLamStar && fLam2 > fLamStar) {

			if (fLam1 >= fLam2) {
				return 0;
			} else {
				return 1;
			}
		}
		System.err.println("argMaxIndex = -1 " + lam1 + "   " + lam2 + "   " + lamStar);
		return -1;
	}

	public double[] calcLineSubOptimal(double[] lineStar, double[] lineMin, double[] lineMax) {

		double[] lineSubOptimal = new double[3];

		if (lineStar[0] >= lineMax[0]) {
			return (lineChars(lineMax[0], 3));
		}

		else if (lineStar[0] <= lineMin[0]) {
			return (lineChars(lineMin[0], 2));
		}

		else {

			double[] lambdaStar12 = calcLambdaStar12(lineStar, lineMin, lineMax);
			double lambdaTemp = argMax(lineStar[0], lambdaStar12[0], lambdaStar12[1]);

			/* lambdaTemp = lambdaStar */
			if (isEqual(lambdaTemp, lineStar[0], precisionLambda)) {
				return (lineStar);
			}

			/*
			 * lambdaStar is NOT subOptimal, so lambdaStar-new will be lambdaTemp and
			 * lambdaMinNew and lambdaMaxNew will be lambdaStar1 or lambdaStar2
			 */
			else {
				/* lambdaTemp == lambdaStar*2 */
				if (isEqual(lambdaTemp, lambdaStar12[1], precisionLambda)) {

					/* Update the Minimum of the interval for the next iteration */
					double[] lineMinNew = lineStar;

					double[] lineStar2 = lineChars(lambdaStar12[1], 2);
					double[] lineTemp = lineStar2;

					lineSubOptimal = calcLineSubOptimal(lineTemp, lineMinNew, lineMax);
				}

				if (isEqual(lambdaTemp, lambdaStar12[0], precisionLambda)) {

					/* Update the Maximum of the interval for the next iteration */
					double[] lineMaxNew = lineStar;

					double[] lineStar1 = lineChars(lambdaStar12[0], 2);
					double[] lineTemp = lineStar1;

					lineSubOptimal = calcLineSubOptimal(lineTemp, lineMin, lineMaxNew);
				}
			}
		}

//		System.out.println("end of calcLambdaSubOpt   "+ lineSubOptimal[0]); 
		return lineSubOptimal;
	}

	public double localMoves(double lambdaSubOptimal, int fb, double dLocalMoves) {

//		System.out.println("begining of local moves   "+lambdaSubOptimal);

		if (lambdaSubOptimal <= minInterval) {
			return minInterval;
		} else if (lambdaSubOptimal >= maxInterval) {
			return maxInterval;
		} else {

			double[] lineSubOptimal = lineChars(lambdaSubOptimal, fb);

			if (isEqual(lineSubOptimal[2], 0, toleranceDerivative)) {
				return lambdaSubOptimal;
			}

			double h = dLocalMoves;
			if (lineSubOptimal[2] < 0) {
				h = (-1) * dLocalMoves;
			}

			double lambdaNew = bringInRange(lambdaSubOptimal + h, minInterval, maxInterval);

			final double fLambdaNew = f(lambdaNew);
			final double fLambdaSubOptimal = f(lambdaSubOptimal);

			if (isEqual(fLambdaNew, fLambdaSubOptimal, toleranceF)) {
				return lambdaSubOptimal;
			}

			else if ((fLambdaNew - fLambdaSubOptimal) > toleranceF) {
				return localMoves(lambdaNew, fb, dLocalMoves);
			}

			else if ((fLambdaSubOptimal - fLambdaNew) > toleranceF) {
				return localMoves(lambdaSubOptimal, fb, dLocalMoves * 0.5);
			}

			/* else do nothing = return the same lambda_subOptimal of last iteration */
			return lambdaSubOptimal;
		}
	}

	public void findOptimalLambda() {
		double[] lineMin = lineChars(minInterval, 2);
		double[] lineMax = lineChars(maxInterval, 3);

		double lambdaStar = findInterSectionLines(lineMin, lineMax);
		double[] lineStar = lineChars(lambdaStar, 2);

		lineStar = lineChars(lambdaStar, 2);

//		System.out.println("inFindOptimalLambda   "+lineStar[0]+"\t"+lineMin[0]+"\t"+lineMax[0]);

		double[] lineSubOptimal = calcLineSubOptimal(lineStar, lineMin, lineMax);

//		System.out.println("inFindOptimalLambda, subOptimalLambda   "+lineSubOptimal[0]);

		this.optimalLambda = localMoves(lineSubOptimal[0], 2, deltaLocalMoves);

//		System.out.println("in findOptimalLamda, optimal lambda: "+optimalLambda);
		this.fOptimalLambda = this.f(optimalLambda);

	}

	public double bringInRange(double lam, double minInt, double maxInt) {
		if (lam <= minInt) {lam = minInt;} 
		else if (lam >= maxInt) {lam = maxInt;}
		return lam;
	}

	public double[] solveNaN(double lam, double minInt, double maxInt) {

		double[] resolvedNan = new double[3];
		if (Double.isNaN(lam) && !Double.isNaN(minInt) && !Double.isNaN(maxInt)) {
			lam = (minInt + maxInt) / 2;
		}
		if (Double.isNaN(minInt) && !Double.isNaN(lam) && !Double.isNaN(maxInt)) {
//			minInt = minInterval + 0.1;
			minInt = lam - 1;
		}
		if (Double.isNaN(maxInt) && !Double.isNaN(minInt) && !Double.isNaN(lam)) {
			maxInt = lam + 1;
		}
		if (Double.isNaN(maxInt) && Double.isNaN(minInt) && Double.isNaN(lam)) {
			System.err.println("ERROR!  " + maxInt + "\t" + minInt + "\t" + lam);
		}

		resolvedNan[0] = bringInRange(lam, 0.01, 9.99);
		resolvedNan[1] = bringInRange(minInt, 0.01, 9.99);
		resolvedNan[2] = bringInRange(maxInt, 0.01, 0.99);

		return resolvedNan;
	}

//	public static List<String> speNamesToDelete(List<String> allSpeNames, int numOfSpeToDelete, int totalNumSpe) {
//		Random rand = new Random();
//		List<String> spe_del = new ArrayList<String>();
//
//		for (int i = 0; i < numOfSpeToDelete; i++) {
//
//			/* random = from + rndGenerator.nextInt(to - from + 1) */
//			int r = rand.nextInt(totalNumSpe);
//
//			if (!spe_del.contains(allSpeNames.get(r))) {
//				spe_del.add(allSpeNames.get(r));
//			} else {
//				numOfSpeToDelete += 1;
//			}
//		}
//		return spe_del;
//	}

	// public static void main(String[] args) {
	//
	// // long startTime = System.currentTimeMillis();
	//// double stepSize = 1e-4;// step calculating derivative
	// double stepSize = 1e-3;
	//// double deltaLocalMoves = 1e-1;
	// double deltaLocalMoves = 1e-2;
	// double tolD = 1e-3;
	// double tolF = 1e-4;
	// double minInterval = 1e-2;
	// double maxInterval = 9.9999;
	// double minAllowed = 1e-2;
	// double maxAllowed = 9.9999;
	// double precisionLambda = 1e-5; // one digit more than the number of digits we
	// require accurately
	// double partitionSize = 0.1;
	// int defaultmaxNodeSize = 100;
	//
	// // int ignoreLine= Integer.parseInt(args[4]);
	//
	// Node root = SpeciesTreeParser.buildInsertWGDsandPartitionTree(args[0],
	// args[1], partitionSize,
	// defaultmaxNodeSize);
	//
	// // Node root = SpeciesTreeParser.buildAndPartitionTree_ReverseEng(args[0],
	// // args[1],
	// // partitionSize, defaultmaxNodeSize,ignoreLine);
	//
	// // args: 0= tree, 1= wgd file, 2= gf count file, 3= number of gf, 4=wgm
	// // retention rates
	// // Node root =
	// SpeciesTreeParser.buildAndPartitionTree_customMultFactor(args[0],
	// // args[1], args[4],
	// // partitionSize, defaultmaxNodeSize);
	// ArrayList<Node> leaves = root.getLeaves();
	//
	// ReadGFcountsFile rgf = new ReadGFcountsFile();
	// List<List<Integer>> gfCounts = rgf.read_all(args[2]);
	//
	// ArrayList<String> gfIDs = rgf.getGfIDs();
	//
	// /** For the 10 Eudicots pruned out of 14 Angio */
	// // double wgtBlen = 0.260;
	// // double oldrootBlen = 0.275776;
	//
	// /** For the 28 Eudicots pruned out of 37 Angio */
	// // double wgtBlen = 0.155383;
	// // double oldrootBlen = 0.2757765996999999;
	//
	// /** For the 8 Monocots pruned out of 37 Angio */
	// // double wgtBlen = 0.222; // it's not really WGT but WGD , just for
	// programming
	// // names
	// // double oldrootBlen = 0.323803;
	//
	// // String idealProfilePath =
	// //
	// "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ComparisonsWithIdealprofiles/IdealProfiles_R1";
	// // ArrayList<Integer> idealProfileR1 =
	// cmmFunct.readColX_Int(idealProfilePath,
	// // 4); // Tau-Mon2_Musa3
	//
	// // int numberOfObservations = 1000;
	// int gf = Integer.parseInt(args[3]);
	// System.out.println(gfIDs.get(Integer.parseInt(args[3])));
	//
	// /** Delete a few species */
	//// String pathToSpeNames=args[4];
	// CommonFunctions cmf = new CommonFunctions();
	//// List<String> allSpeciesNames = cmf.read1ColFile_String(pathToSpeNames);
	//
	//// int numOfSpeToDelete=Integer.parseInt(args[5]);
	//// int numOfSpeToDelete = 5;
	//// List<String> speNamesToDelete = speNamesToDelete(allSpeciesNames,
	// numOfSpeToDelete, allSpeciesNames.size());
	//
	//
	// for (int rootSize = 1; rootSize <= 10; rootSize++) {
	// // int rootSize= Integer.parseInt(args[3]);
	//
	// /**
	// * Don't move these lines out of the loop. because pvalues change the leaf
	// * values and you have to reset once you move to the next root size within the
	// * same gene family: for now we don't calculate p-values
	// */
	//
	// SpeciesTreeParser.setLeavesValues(root, gfCounts, gf);
	// // SpeciesTreeParser.setLeavesValues_one(root, idealProfileR1);
	//// root.deleteNodeFromTree(speNamesToDelete);
	// ArrayList<Node> speciesTree = SpeciesTreeParser.setMaxNodeSize(root,
	// defaultmaxNodeSize);
	//
	// ProbCalculator probCache = new ProbCalculator();
	//
	// CuttingPlaneMethod cpm = new CuttingPlaneMethod(speciesTree, rootSize,
	// stepSize, deltaLocalMoves, tolD,
	// tolF, minInterval, maxInterval, minAllowed, maxAllowed, precisionLambda,
	// probCache, root, gfCounts);
	//
	// // RootSizeThread runThread = new RootSizeThread(cpm);
	// // runThread.start();
	//
	// cpm.findOptimalLambda();
	//
	// System.out.print(rootSize + "\t" + cpm.getOptimalLambda() + "\t" +
	// cpm.getfOptimalLambda());
	//
	// // Pvalues pv = new Pvalues(root, speciesTree, cpm.getOptimalLambda(),
	// // numberOfObservations, probCache);
	//
	// // int eudicotSize =
	// // IncludeAngioSperms.generateSizeAtEudicots(root.getmaxNodeSize(),
	// // rootSize, cpm.getOptimalLambda(), wgtBlen, oldrootBlen,
	// // probCache);
	// // double pValue = pv.calculateConditionalPvalues(rootSize,
	// // eudicotSize, cpm.getfOptimalLambda());
	//
	// // double pValue = pv.calculateConditionalPvalues(rootSize, rootSize,
	// // cpm.getfOptimalLambda());
	// //
	// // System.out.print("\t" + pValue);
	//
	// System.out.print("\n");
	//
	// }
	// // long endTime = System.currentTimeMillis();
	// // System.out.println(endTime - startTime);

	// }

}
