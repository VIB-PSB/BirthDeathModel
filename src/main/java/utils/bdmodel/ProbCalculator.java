package utils.bdmodel;

import java.util.Map;
import java.util.TreeMap;

public class ProbCalculator {

	/**
	 * Key for the caching {@link Map} that maps the four parameters lambda,
	 * branch length, parent size and child size onto its probability.
	 *
	 */
	public class CacheKey implements Comparable<CacheKey> {

		static public final double PRECISION = 0.00001;

		final double t;
		final int s;
		final int c;

		public CacheKey(double t2, int s2, int c2) {
			t = t2;
			s = s2;
			c = c2;
		}

		private boolean precEqual(double x, double y) {
			return Math.abs(x - y) < PRECISION;
		}

		@Override
		public int compareTo(CacheKey rhs) {
			if (!precEqual(t, rhs.t))
				return t > rhs.t ? 1 : -1;
			if (s != rhs.s)
				return s > rhs.s ? 1 : -1;
			if (c != rhs.c)
				return c > rhs.c ? 1 : -1;
			return 0;
		}

		@Override
		public boolean equals(Object obj) {
			if (!(obj instanceof CacheKey))
				return false;
			CacheKey rhs = (CacheKey) obj;
			return this.compareTo(rhs) == 0;
		}
	}

	private static final boolean CACHING_ENABLED = true;

	private Map<CacheKey, Double> cache;
	private double lastLambda;
	private int cacheHit;
	private int cacheMiss;

	public ProbCalculator() {
		cache = new TreeMap<CacheKey, Double>();
	}

	/**
	 * Create a new ProbCalculator using a cache from a former iteration
	 * 
	 * @param cache
	 */
	public ProbCalculator(TreeMap<CacheKey, Double> cache) {
		this.cache = cache;
		cacheHit = 0;
		cacheMiss = 0;
	}

	public void setCache(TreeMap<CacheKey, Double> cache) {
		this.cache = cache;
	}

	/**
	 * @param lambda
	 *            : gene duplication rate
	 * @param t
	 *            : branch length
	 * @param s
	 *            : parent size
	 * @param c
	 *            : child size
	 * @return transition probability of going from s gene counts at a parent
	 *         node to c gene counts for the child node, given Bailey's BD
	 *         model, 1964.
	 */
	public double probCalc(double lambda, double t, int s, int c) {
		if (CACHING_ENABLED) {
			if (cache != null) {
				return probCalcCache(lambda, t, s, c);
			}
			// else {
			// System.err.println("WARNING: NOT USING CACHING!");
			// }
		}

		if (lambda < 0)
			System.out.println("Error: negative lambda passed: " + lambda);
		if (t < 0)
			System.out.println("Error: negative branch length passed: " + t);
		if (s < 0)
			System.out.println("Error: invalid parent size passed: " + s);
		if (c < 0)
			System.out.println("Error: invalid child size passed: " + c);

		int b = Math.min(s, c);

		double alpha = (lambda * t) / (1 + lambda * t);

		double sum = 0;

		for (int j = 0; j < b + 1; j++) {

			double a = MathOperations.binomial(s, j);

			double d = MathOperations.binomial((s + c - j - 1), (s - 1));

			double e = Math.pow(alpha, (s + c - (2 * j)));

			double h = Math.pow((1 - (2 * alpha)), j);

			double res = a * d * e * h;
			sum += res;
		}

		return sum;
	}

	private double probCalcCache(double lambda, double t, int s, int c) {
		if (lambda != lastLambda) {
			cache.clear();
			lastLambda = lambda;
		}
		CacheKey key = new CacheKey(t, s, c);
		Double cacheVal = cache.get(key);
		if (cacheVal != null) {
			cacheHit++;
			return cacheVal;
		}
		cacheMiss++;

		int b = Math.min(s, c);

		double alpha = (lambda * t) / (1 + lambda * t);

		double sum = 0;

		for (int j = 0; j < b + 1; j++) {

			double a = MathOperations.binomial(s, j);

			double d = MathOperations.binomial(s + c - j - 1, s - 1);

			double e = Math.pow(alpha, (s + c - (2 * j)));

			double h = Math.pow((1 - (2 * alpha)), j);

			double res = a * d * e * h;
			sum += res;
		}
		cache.put(key, sum);
		return sum;
	}
}
