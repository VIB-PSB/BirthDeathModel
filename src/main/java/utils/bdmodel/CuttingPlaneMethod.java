package utils.bdmodel;


import java.util.HashMap;
import java.util.List;


public class CuttingPlaneMethod {

	private List<Node> aryln;
	private int rootSize;
	private Node root;

	private final double stepSize;
	private final double deltaLocalMoves;
	private final double toleranceDerivative;
	private final double toleranceF;
	private final double precisionLambda;
	private final double minInterval;
	private final double maxInterval;
	public int counter;
	private double optimalLambda;
	private double fOptimalLambda;
	private ProbCalculator probCalc;
	private HashMap<Double, Double> fCache;

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
			double tolD, double tolF, double minInterval, double maxInterval,
			double precision, ProbCalculator probCache, Node root) {

		this.stepSize = stepSize;
		this.deltaLocalMoves = deltaLocalMoves;
		this.rootSize = rootSize;
		this.toleranceDerivative = tolD;
		this.toleranceF = tolF;
		this.minInterval = minInterval;
		this.maxInterval = maxInterval;
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

	public double forwardDiff(double lambda) {

		double fLambda = f(lambda);
		double fLambdaPlusDelta = f(lambda + stepSize);

		double forwardDiff = (fLambdaPlusDelta - fLambda) / stepSize;

		return forwardDiff;
	}

	public double backwardDiff(double lambda) {

		double fLambda = f(lambda);
		double fLambdaMinusDelta = f(lambda - stepSize);

		double backwardDiff = (fLambda - fLambdaMinusDelta) / stepSize;

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
		double lambdaStar = ((line1[2] * line1[0]) - (line2[2] * line2[0]) + line2[1] - line1[1])
				/ (line1[2] - line2[2]);

		double[] noNaN = solveNaN(lambdaStar, line1[0], line2[0]);
		return noNaN[0]; // lambdaStar_edited
	}

	/* returns [lam*_1 , lam*_2] */
	public double[] calcLambdaStar12(double[] lineStar, double[] line1, double[] line2) {
		double[] noNaN = solveNaN(lineStar[0], line1[0], line2[0]);

		lineStar = lineChars(noNaN[0], 2);
		line1 = lineChars(noNaN[1], 2);
		line2 = lineChars(noNaN[2], 2);

		double[] lambdaStar12 = new double[2];
		double lambdaStar_1 = findInterSectionLines(lineStar, line1);
		double lambdaStar_2 = findInterSectionLines(lineStar, line2);

		lambdaStar12[0] = lambdaStar_1;
		lambdaStar12[1] = lambdaStar_2;

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
			return lamStar; // Sub-optimality
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
			return 2; // Sub-optimality
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

			if (isEqual(lambdaTemp, lineStar[0], precisionLambda)) {
				return (lineStar);
			}

			else {
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

		return lineSubOptimal;
	}

	public double localMoves(double lambdaSubOptimal, int fb, double dLocalMoves) {

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

		double[] lineSubOptimal = calcLineSubOptimal(lineStar, lineMin, lineMax);
		this.optimalLambda = localMoves(lineSubOptimal[0], 2, deltaLocalMoves);
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
}
