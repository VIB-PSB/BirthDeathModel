package utils.bdmodel;

import cern.jet.stat.Gamma;

public class MathOperations {

	static public final boolean CACHING_ENABLED = true;
	static public final int CACHE_SIZE = 500;
	
	static public double cache[][] = generateCache(CACHE_SIZE);

	static private double[][] generateCache(int N) {
		double mycache[][] = new double[N][N];
		mycache[0][0] = 1;
		for (int n = 1; n < N; ++n) {
			for (int k = 0; k < n + 1; ++k) {
				mycache[n][k] = binomialCalc(n, k);
			}
		}
		return mycache;
	}

	static public double binomial(int n, int k) {
		if (CACHING_ENABLED && n<CACHE_SIZE) {
			return cache[n][k];
		} else {
			return binomialCalc(n, k);
		}
	}

	/**
	 * @return: binomial coefficient k out of n
	 */

	static private double binomialCalc(int n, int k) {

		if ((n==0 && k!=n) || k>n){
			throw new RuntimeException("Weird binomials!");
		}
		
		if(k==0 || k==n){
			return 1;
		}
		
		else if (k==1 || k== n-1 ){
			return n;
		}
		
		else if(k==2 || k== n-2){
			return (n*(n-1)/2);
		}
		
		else {

			double a = Gamma.logGamma(n + 1);

			double b = Gamma.logGamma(k + 1);

			double c = Gamma.logGamma(n - k + 1);

			double d = Math.exp(a - b - c);

			// return (Math.exp(Gamma.logGamma(n + 1) - Gamma.logGamma(k + 1)
			// - Gamma.logGamma(n - k + 1)));

			return d;
		}

	}

	/**
	 * @param args
	 */
	 public static double [] giveLogArray(double [] a){
		
		double [] b= new double[a.length];
		
		for(int i=0; i< a.length; i++){
			
			b[i] = cern.jet.math.Arithmetic.log10(a[i]);
		}
		
		return b;
	}	 
}
