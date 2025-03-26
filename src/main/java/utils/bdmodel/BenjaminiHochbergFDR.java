package utils.bdmodel;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;

/**
 * * Copyright (c) 2005 Flanders Interuniversitary Institute for Biotechnology (VIB)
 * *
 * * Authors : Steven Maere, Karel Heymans
 * *
 * * This program is free software; you can redistribute it and/or modify
 * * it under the terms of the GNU General Public License as published by
 * * the Free Software Foundation; either version 2 of the License, or
 * * (at your option) any later version.
 * *
 * * This program is distributed in the hope that it will be useful,
 * * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * * The software and documentation provided hereunder is on an "as is" basis,
 * * and the Flanders Interuniversitary Institute for Biotechnology
 * * has no obligations to provide maintenance, support,
 * * updates, enhancements or modifications.  In no event shall the
 * * Flanders Interuniversitary Institute for Biotechnology
 * * be liable to any party for direct, indirect, special,
 * * incidental or consequential damages, including lost profits, arising
 * * out of the use of this software and its documentation, even if
 * * the Flanders Interuniversitary Institute for Biotechnology
 * * has been advised of the possibility of such damage. See the
 * * GNU General Public License for more details.
 * *
 * * You should have received a copy of the GNU General Public License
 * * along with this program; if not, write to the Free Software
 * * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * *
 * * Authors: Steven Maere, Karel Heymans
 * * Date: Mar.25.2005
 * * Description: Class implementing the Benjamini and Hochberg FDR correction algorithm.        
 **/

/**
 * ************************************************************************
 * BenjaminiHochbergFDR.java: Steven Maere & Karel Heymans (c) March 2005
 * ------------------------
 * <p/>
 * Class implementing the Benjamini and Hochberg FDR correction algorithm.
 * <p/>
 * ************************************************************************
 */


// Class Declaration
public class BenjaminiHochbergFDR {

	// Fields (instance variables): attributes or properties of the class
	/**
	 * the raw p-values that were given as input for the constructor, order
	 * corresponds to String [] goLabels.
	 */
	private double[] pvalues;
	/**
	 * the test IDs ordened according to the ordened pvalues.
	 */
	private int[] ordenedIndices;
	/**
	 * the raw p-values ordened in ascending order.
	 */
	private double[] ordenedPvalues;
	/**
	 * the adjusted p-values ordened in ascending order.
	 */
	private double[] adjustedOrdenedPvalues;

	/**
	 * results (adjusted p-values) in original order.
	 */
	private double[] adjustedPvalues;

	/**
	 * the significance level.
	 */
	private BigDecimal alpha;
	/**
	 * the number of tests.
	 */
	private int m;
	/**
	 * scale for the division in de method 'runFDR'.
	 */
	private static final int RESULT_SCALE = 100;


	/*--------------------------------------------------------------
	CONSTRUCTOR.
	--------------------------------------------------------------*/

	/**
	 * Constructor.
	 * 
	 * @param pvalues
	 *            Hashmap of Strings with the test Labels and their pvalues.
	 * @param alpha
	 *            String with the desired significance level.
	 */

	public BenjaminiHochbergFDR(double[] pvalues, double alpha) {
		// Get all the go labels and their corresponding pvalues from the map

		this.pvalues = pvalues;
		this.alpha = new BigDecimal(alpha);
		this.m = pvalues.length;
        this.ordenedPvalues = new double[m];
        this.ordenedIndices = new int[m];
		this.adjustedOrdenedPvalues = new double[m];
        this.adjustedPvalues = new double[m];
	}

	/*--------------------------------------------------------------
	METHODS.
	--------------------------------------------------------------*/

	/**
	 * Method that calculates the Benjamini and Hochberg correction of the false
	 * discovery rate NOTE : convert array indexes [0..m-1] to ranks [1..m].
	 * orden raw p-values low .. high test p<(i/m)*alpha from high to low (for
	 * i=m..1) i* (istar) first i such that the inequality is correct. reject
	 * hypothesis for i=1..i* : labels 1..i* are overrepresented
	 * <p/>
	 * adjusted p-value for i-th ranked p-value p_i^adj = min(k=i..m)[min(1,m/k
	 * p_k)]
	 */

  public class Pair implements Comparable<Pair> {
      public final int index;
      public final double value;

      public Pair(int index, double value) {
          this.index = index;
          this.value = value;
      }

      @Override
      public int compareTo(Pair other) {
          return Double.valueOf(this.value).compareTo(other.value);
      }
  }
        
	public void calculate() {
	
		// Ordening the pvalues.
        Pair[] ordenedPvals = new Pair[pvalues.length];
		for (int i=0; i<pvalues.length; i++){
                    ordenedPvals[i]= new Pair(i,pvalues[i]);
                }
        Arrays.sort(ordenedPvals);
        for (int i=0; i<ordenedPvals.length; i++){
                    ordenedPvalues[i] = ordenedPvals[i].value;
                    ordenedIndices[i] = ordenedPvals[i].index;
                }
                
		// Calculating adjusted p-values.
		BigDecimal min = new BigDecimal("" + 1);
		BigDecimal mkprk;
		for (int i = m; i > 0; i--) {
			mkprk = (new BigDecimal("" + m).multiply(new BigDecimal(ordenedPvalues[i - 1]))).divide(new BigDecimal(""
					+ i), RESULT_SCALE, RoundingMode.HALF_UP);
			if (mkprk.compareTo(min) < 0) {
				min = mkprk;
			}
			adjustedOrdenedPvalues[i - 1] = min.doubleValue();

		}
                
		for (int i = 0; i < ordenedIndices.length; i++) {
			adjustedPvalues[ordenedIndices[i]] = adjustedOrdenedPvalues[i];
		}
	}


	/*--------------------------------------------------------------
	  GETTERS.
	--------------------------------------------------------------*/


	public double[] getAdjustedPvalues() {
		return adjustedPvalues;
	}

}
