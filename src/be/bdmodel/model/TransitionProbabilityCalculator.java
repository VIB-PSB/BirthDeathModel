package be.ugent.psb.setas.bdmodel.model;

/*
 * #%L
 * BirthDeathModel
 * %%
 * Copyright (C) 2017 VIB/PSB/UGent - Setareh Tasdighian
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

import java.util.Map;
import java.util.TreeMap;

public class TransitionProbabilityCalculator {

	/**
	 * Key for the caching {@link Map} that maps the four parameters lambda,
	 * branch length (t), parent size (s) and child size (c) onto its probability.
	 *
	 */
	public class CacheKey implements Comparable<CacheKey> {

		static public final double PRECISION = 1e-6;

		final double t;
		final int s;
		final int c;

		public CacheKey(double t2, int s2, int c2) {
			t = t2;
			s = s2;
			c = c2;
		}

		private boolean precisionEquality(double x, double y) {
			return Math.abs(x - y) < PRECISION;
		}

		@Override
		public int compareTo(CacheKey rhs) {
			if (!precisionEquality(t, rhs.t))
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

	public TransitionProbabilityCalculator() {
		cache = new TreeMap<CacheKey, Double>();
	}

	/**
	 * Create a new TransitionProbabilityCalculator using a cache from a former iteration
	 * 
	 * @param cache
	 */
	public TransitionProbabilityCalculator(TreeMap<CacheKey, Double> cache) {
		this.cache = cache;
	}

	public void setCache(TreeMap<CacheKey, Double> cache) {
		this.cache = cache;
	}

	/**
	 * @param lambda
	 *            : gene duplication rate
	 * @param branchLength
	 *            : branch length
	 * @param parentGeneCount
	 *            : parent size
	 * @param childGeneCount
	 *            : child size
	 * @return transition probability of going from s gene counts at a parent
	 *         node to c gene counts for the child node, given baily's BD
	 *         be.ugent.psb.setas.bdmodel.model, 1964.
	 */
	public double probabilityCalculator(double lambda, double branchLength, int parentGeneCount, int childGeneCount) {
		if (CACHING_ENABLED) {
			if (cache != null) {
				return probabilityCalculatorCache(lambda, branchLength, parentGeneCount, childGeneCount);
			}
		}

		if (lambda < 0)
			System.err.println("Error: negative value of lambda is passed: " + lambda);
		if (branchLength < 0)
			System.err.println("Error: negative value for branch length is passed: " + branchLength);
		if (parentGeneCount < 0)
			System.err.println("Error: negative parent size is passed: " + parentGeneCount);
		if (childGeneCount < 0)
			System.err.println("Error: negative child size is passed: " + childGeneCount);

		int minSandC = Math.min(parentGeneCount, childGeneCount);

		double alpha = (lambda * branchLength) / (1 + lambda * branchLength);

		double sum = 0;

		for (int j = 0; j < minSandC + 1; j++) {

			double a = MathematicalOperations.binomial(parentGeneCount, j);
			double d = MathematicalOperations.binomial((parentGeneCount + childGeneCount - j - 1), (parentGeneCount - 1));
			double e = Math.pow(alpha, (parentGeneCount + childGeneCount - (2 * j)));
			double h = Math.pow((1 - (2 * alpha)), j);

			double result = a * d * e * h;
			sum += result;
		}

		return sum;
	}

	private double probabilityCalculatorCache(double lambda, double branchLength, int parentGeneCount, int childGeneCount) {
		if (lambda != lastLambda) {
			cache.clear();
			lastLambda = lambda;
		}
		CacheKey key = new CacheKey(branchLength, parentGeneCount, childGeneCount);
		Double cacheVal = cache.get(key);
		if (cacheVal != null) {
			return cacheVal;
		}
		int minSandC = Math.min(parentGeneCount, childGeneCount);

		double alpha = (lambda * branchLength) / (1 + lambda * branchLength);

		double sum = 0;

		for (int j = 0; j < minSandC + 1; j++) {

			double a = MathematicalOperations.binomial(parentGeneCount, j);

			double d = MathematicalOperations.binomial(parentGeneCount + childGeneCount - j - 1, parentGeneCount - 1);

			double e = Math.pow(alpha, (parentGeneCount + childGeneCount - (2 * j)));

			double h = Math.pow((1 - (2 * alpha)), j);

			double result = a * d * e * h;
			sum += result;
		}
		cache.put(key, sum);
		return sum;
	}

}
