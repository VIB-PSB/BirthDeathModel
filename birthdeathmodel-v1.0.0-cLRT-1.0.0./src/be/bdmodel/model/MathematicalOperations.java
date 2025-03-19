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

import cern.jet.stat.Gamma;

/**
 * Using the library Cern, this class helps calcuate the binomial coefficients,
 * and the log(10) of a given array of doubles.
 * 
 * @author setas
 *
 */
public class MathematicalOperations {

	static public final boolean CACHING_ENABLED = true;
	static public final int CACHE_SIZE = 500;

	static public double cache[][] = generateCache(CACHE_SIZE);

	static private double[][] generateCache(int my_integer) {
		double mycache[][] = new double[my_integer][my_integer];
		mycache[0][0] = 1;
		for (int n = 1; n < my_integer; ++n) {
			for (int k = 0; k < n + 1; ++k) {
				mycache[n][k] = binomialCalc(n, k);
			}
		}
		return mycache;
	}

	static public double binomial(int n, int k) {
		if (CACHING_ENABLED && n < CACHE_SIZE) {
			return cache[n][k];
		} else {
			return binomialCalc(n, k);
		}
	}

	/**
	 * @return: binomial coefficient k out of n
	 */
	static private double binomialCalc(int n, int k) {

		if ((n == 0 && k != n) || k > n) {
			throw new RuntimeException("Calculating binomials for n = "+n+" k= "+k);
		}

		if (k == 0 || k == n) {
			return 1;
		}

		else if (k == 1 || k == n - 1) {
			return n;
		}

		else if (k == 2 || k == n - 2) {
			return (n * (n - 1) / 2);
		}

		else {

			double a = Gamma.logGamma(n + 1);

			double b = Gamma.logGamma(k + 1);

			double c = Gamma.logGamma(n - k + 1);

			double d = Math.exp(a - b - c);

			return d;
		}

	}
	
	/**
	 * Calculating log(10) of values in an array
	 */
	public static double[] giveLogarithm10Array(double[] arrayOfdoubles) {

		double[] logarithm10Array = new double[arrayOfdoubles.length];

		for (int i = 0; i < arrayOfdoubles.length; i++) {

			logarithm10Array[i] = cern.jet.math.Arithmetic.log10(arrayOfdoubles[i]);
		}

		return logarithm10Array;
	}

}
