package workflows.negativeWGMs;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import utils.bdmodel.BenjaminiHochbergFDR;
import utils.parsers.CommonFunctions;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

// Class Declaration
public class NegativeWGM_PostprocessingGFListFDR {

	// Fields (instance variables): attributes or properties of the class
	private final String negatives_H0Full_outputDir;
	private final String negatives_H0Neg_outputDir;
	private final String origLRTDir;
	private final String listGFs; // Comma-separated list of GF numbers to be investigated
	private final int wgdnumber;
	private final double fdr;

	// Constructor to initialize an object of this class
	public NegativeWGM_PostprocessingGFListFDR(String negatives_H0Full_outputDir, String negatives_H0Neg_outputDir, String origLRTDir, String listGFs, int wgdnumber, double fdr) {
		this.negatives_H0Full_outputDir = negatives_H0Full_outputDir;
		this.negatives_H0Neg_outputDir = negatives_H0Neg_outputDir;
		this.origLRTDir = origLRTDir;
		this.listGFs = listGFs;
		this.wgdnumber = wgdnumber;
		this.fdr = fdr;
	}

	// Methods: behaviors or operations that objects of the class can perform
	public static double calculateLRT(double loglkH0, double loglkH1) {
		return (2 * (loglkH1 - loglkH0));
	}

	public void execute() {

		/* ****************************************************************** */
		/* Decision-making rules for complementary LRTs on all negative WGMs
		/* Accepts as input a LIST of GF ranks.
		/* ****************************************************************** */

		DecimalFormat df2 = new DecimalFormat("0.00000");
		CommonFunctions cmf = new CommonFunctions();
		String dirName = System.getProperty("user.dir");
		
		File results_perGF = new File(dirName, "RvsEng_Negatives_perGF_FDR.txt");
		File results_perWGD_Top = new File(dirName, "RvsEng_Negatives_perWGD_Top_FDR.txt");
		results_perGF.delete();
		results_perWGD_Top.delete();
		BufferedWriter output_perGF;
		BufferedWriter output_perWGD_Top;
		
		try {
			results_perGF.createNewFile();
			results_perWGD_Top.createNewFile();
			output_perGF = new BufferedWriter(new FileWriter(results_perGF,true));
			output_perWGD_Top = new BufferedWriter(new FileWriter(results_perWGD_Top,true));

			System.out.println("GF rank (starts from 0)\tGF ID\tnegative WGD number (starts from 1)\toriginal LRT under H0 = negative WGD absent\tmedian of simulated LRT under H0 = negative WGD absent\t95th percentile of simulated LRT under H0 = negative WGD absent\tFDR-corrected pval H0 = negative WGD absent\toriginal LRT under H0 = negative WGD present\tmedian of simulated LRT under H0 = negative WGD present\t95th percentile of simulated LRT under H0 = negative WGD present\tFDR-corrected pval H0 = negative WGD present\tinference");
			
			String gfID = new String();
			String[] GFstring = listGFs.split(",");
			
			int[] GFs = new int[GFstring.length];
			for (int i = 0; i < GFstring.length; i++) {
				GFs[i] = Integer.parseInt(GFstring[i]);
			}

			int[] nrNegWGDsRejected = new int[GFs.length];
			int[] nrNegWGDsAccepted = new int[GFs.length];
			int[] nrUndecided = new int[GFs.length];
			int[] nrBothRejected = new int[GFs.length];
			
			int[] nrTopGFsReject = new int[wgdnumber];
			int[] nrTopGFsAccept = new int[wgdnumber];
			int[] nrTopGFsUndecided = new int[wgdnumber];
			int[] nrTopGFsRejectBoth = new int[wgdnumber];

			for (int negWgd = 1; negWgd<=wgdnumber; negWgd++) {
				nrTopGFsReject[negWgd-1] = 0;
				nrTopGFsAccept[negWgd-1] = 0;
				nrTopGFsUndecided[negWgd-1] = 0;
				nrTopGFsRejectBoth[negWgd-1] = 0;
			}
			
			double[][] pH0Full = new double[GFstring.length][wgdnumber];
			// Bonferroni-corrected pH0Full
			double[][] BonfpH0Full = new double[GFstring.length][wgdnumber];
			// Benjamini-Hochberg corrected pH0Full
			double[][] FDRpH0Full = new double[GFstring.length][wgdnumber];

			double[][] pH0Neg = new double[GFstring.length][wgdnumber];
			// Bonferroni-corrected pH0Rm
			double[][] BonfpH0Neg = new double[GFstring.length][wgdnumber];
			// Benjamini-Hochberg corrected pH0Rm
			double[][] FDRpH0Neg = new double[GFstring.length][wgdnumber];

			for (int k=0; k<GFs.length ; k++) {
				// i = GF number
				int i = GFs[k];
				// Read in LRT test results for real data under H0Full
				// Read in file i from origLRTDir
				String fileName = String.format("output_%d.txt",i);
				String origFile = new File(origLRTDir,fileName).toString();
				List<List<String>> outputOrigLRT = cmf.readMapFile(origFile);
				// GF ID first on first line
				gfID = outputOrigLRT.get(0).get(0);
				// Likelihood under H0 = full tree also on first line, 3 places away from GF ID
				double origLogLkFullTree = Double.parseDouble(outputOrigLRT.get(0).get(3));
				double[] origLogLkNegTree = new double[wgdnumber];
				double[] origLRTH0Full = new double[wgdnumber];
				double[] origLRTH0Neg = new double[wgdnumber];
				for (int negWgd = 1; negWgd<=wgdnumber; negWgd++) {
					origLogLkNegTree[negWgd-1] = Double.parseDouble(outputOrigLRT.get(negWgd).get(2));
					origLRTH0Full[negWgd-1] = calculateLRT(origLogLkFullTree, origLogLkNegTree[negWgd-1]);
					origLRTH0Neg[negWgd-1] = calculateLRT(origLogLkNegTree[negWgd-1], origLogLkFullTree);
				}
				// Read in LRT test results for 1000 simulated observations under H0Full
				double[][] simLRTH0Full = new double[wgdnumber][1000];
				for (int j = 0; j<1000; j++) {
					// Read from 1000 files in negatives_H0Full_outputDir for GF i
					// !! numbering of output files starts at 1 instead of 0
					int filenr = j+1;
					String fileName0 = String.format("Output_Neg_SimH0Full_%d/output_%d.txt",i,filenr);
					String simFile0 = new File(negatives_H0Full_outputDir,fileName0).toString();
					List<List<String>> outputH0Full = cmf.readMapFile(simFile0);
					// First line is GF ID
					String gfID2 = outputH0Full.get(0).get(0);
					if (!gfID.equalsIgnoreCase(gfID2)) {
						System.out.println("error with gfID match");
					}

					for (int negWgd = 1; negWgd<=wgdnumber; negWgd++) {
							simLRTH0Full[negWgd-1][j] = Double.parseDouble(outputH0Full.get(negWgd).get(6));
					}
				}
				// Read in LRT test results for 1000 simulated observations under H0Neg
				double[][] simLRTH0Neg = new double[wgdnumber][1000];
				for (int j = 0; j<1000; j++) {
					// Read from 1000 files in negatives_H0Neg_outputDir for GF i
					// !! numbering of output files starts at 1 instead of 0
					int filenr = j+1;
					String fileName1 = String.format("Output_Neg_SimH0Neg_%d/output_%d.txt",i,filenr);
					String simFile1 = new File(negatives_H0Neg_outputDir,fileName1).toString();
					List<List<String>> outputH0Neg = cmf.readMapFile(simFile1);
					// First line is GF ID
					String gfID2 = outputH0Neg.get(0).get(0);
					if (!gfID.equalsIgnoreCase(gfID2)) {
						System.out.println("error with gfID match");
					}

					for (int negWgd = 1; negWgd<=wgdnumber; negWgd++) {
						simLRTH0Neg[negWgd-1][j] = Double.parseDouble(outputH0Neg.get(negWgd).get(6));
					}
				}
				
				double[] medianH0FullDist = new double[wgdnumber];
				double[] medianH0NegDist = new double[wgdnumber];
				double[] pct95H0FullDist = new double[wgdnumber];
				double[] pct95H0NegDist = new double[wgdnumber];
				
				// Calculate percentiles of origLRTH0Full in simLRTH0Full and origLRTH0Neg in simLRTH0Neg
				nrNegWGDsRejected[k]=0;
				nrNegWGDsAccepted[k]=0;
				nrUndecided[k]=0;
				nrBothRejected[k]=0;

				for (int negWgd = 1; negWgd<=wgdnumber; negWgd++) {
					Percentile p = new Percentile();
					Median m = new Median();
					double[] simLRTH0FullDist = simLRTH0Full[negWgd-1];
					// Sort from low to high
					Arrays.sort(simLRTH0FullDist);
					medianH0FullDist[negWgd-1] = m.evaluate(simLRTH0FullDist);
					pct95H0FullDist[negWgd-1] = p.evaluate(simLRTH0FullDist, 95);
					double[] simLRTH0NegDist = simLRTH0Neg[negWgd-1];
					// nH0Full = number of times model on permuted data for H0 = full tree gives equal or better LRT than model H0 = full tree on original data
					double nH0Full = 0;
					// Number of permutations
					double kH0Full = 1000;
					for (int j = 0; j < simLRTH0Full[negWgd-1].length; j++) {
						if (origLRTH0Full[negWgd-1]<=simLRTH0Full[negWgd-1][j]) {
							nH0Full = nH0Full + 1;
						}
					}
					pH0Full[k][negWgd-1] = (nH0Full+1)/(kH0Full+1);
					BonfpH0Full[k][negWgd-1] = pH0Full[k][negWgd-1]*wgdnumber;

					Arrays.sort(simLRTH0NegDist);
					medianH0NegDist[negWgd-1] = m.evaluate(simLRTH0NegDist);
					pct95H0NegDist[negWgd-1] = p.evaluate(simLRTH0NegDist, 95);
					// nH0Neg = number of times model on permuted data for H0 = neg WGD present gives equal or better LRT than model H0 = neg WGD present on original data
					double nH0Rm = 0;
					// Number of permutations
					double kH0Rm = 1000;
					for (int j = 0; j < simLRTH0Neg[negWgd-1].length; j++) {
						if (origLRTH0Neg[negWgd-1]<=simLRTH0Neg[negWgd-1][j]) {
							nH0Rm = nH0Rm + 1 ;
						}
					}
					pH0Neg[k][negWgd-1] = (nH0Rm+1)/(kH0Rm+1);
					BonfpH0Neg[k][negWgd-1] = pH0Neg[k][negWgd-1]*wgdnumber;
				}

				BenjaminiHochbergFDR fdr1 = new BenjaminiHochbergFDR(pH0Full[k],fdr);
				fdr1.calculate();
				FDRpH0Full[k] = fdr1.getAdjustedPvalues();

				BenjaminiHochbergFDR fdr2 = new BenjaminiHochbergFDR(pH0Neg[k],fdr);
				fdr2.calculate();
				FDRpH0Neg[k] = fdr2.getAdjustedPvalues();

				for (int negWgd = 1; negWgd<=wgdnumber; negWgd++) {
					// If actual LRT value higher than 95th percentile (i.e. worse likelihood of H0 versus H1 in reality than in simulations under H0) --> reject H0
					// !! take into account that numbering of GFs starts from 0 for final output, WGD numbering starts at 1
					int gfNumber = i;
					if (FDRpH0Full[k][negWgd-1] < fdr) {
						if (FDRpH0Neg[k][negWgd-1] < fdr) {
							// Both H0 rejected --> neither model is appropriate
							System.out.println(gfNumber + "\t" + gfID + "\t" + negWgd + "\t" + df2.format(origLRTH0Full[negWgd-1]) + "\t" + df2.format(medianH0FullDist[negWgd-1]) + "\t" + df2.format(pct95H0FullDist[negWgd-1]) + "\t" + df2.format(FDRpH0Full[k][negWgd-1]) + "\t" + df2.format(origLRTH0Neg[negWgd-1]) + "\t" + df2.format(medianH0NegDist[negWgd-1]) + "\t" + df2.format(pct95H0NegDist[negWgd-1]) + "\t" + df2.format(FDRpH0Neg[k][negWgd-1]) + "\tboth H0 rejected");
							nrBothRejected[k]++;
							nrTopGFsRejectBoth[negWgd-1]++;
						}
						else {
							// H0 = full tree rejected and H0 = presence of neg WGD not rejected --> presence of negative WGD preferred
							System.out.println(gfNumber + "\t" + gfID + "\t" + negWgd + "\t" + df2.format(origLRTH0Full[negWgd-1]) + "\t" + df2.format(medianH0FullDist[negWgd-1]) + "\t" + df2.format(pct95H0FullDist[negWgd-1]) + "\t" + df2.format(FDRpH0Full[k][negWgd-1]) + "\t" + df2.format(origLRTH0Neg[negWgd-1]) + "\t" + df2.format(medianH0NegDist[negWgd-1]) + "\t" + df2.format(pct95H0NegDist[negWgd-1]) + "\t" + df2.format(FDRpH0Neg[k][negWgd-1]) + "\tnegative WGD present");
							nrNegWGDsAccepted[k]++;
							nrTopGFsAccept[negWgd-1]++;
						}
					}
					else {
						if (FDRpH0Neg[k][negWgd-1] < fdr) {
							// H0 = full tree not rejected and H0 = presence of neg WGD rejected --> absence of negative WGD preferred
							System.out.println(gfNumber + "\t" + gfID + "\t" + negWgd + "\t" + df2.format(origLRTH0Full[negWgd-1]) + "\t" + df2.format(medianH0FullDist[negWgd-1]) + "\t" + df2.format(pct95H0FullDist[negWgd-1]) + "\t" + df2.format(FDRpH0Full[k][negWgd-1]) + "\t" + df2.format(origLRTH0Neg[negWgd-1]) + "\t" + df2.format(medianH0NegDist[negWgd-1]) + "\t" + df2.format(pct95H0NegDist[negWgd-1]) + "\t" + df2.format(FDRpH0Neg[k][negWgd-1]) + "\tnegative WGD absent");
							nrNegWGDsRejected[k]++;
							nrTopGFsReject[negWgd-1]++;
						}
						else {
							// H0 = full tree not rejected and H0 = presence of neg WGD not rejected --> no preference for either model
							System.out.println(gfNumber + "\t" + gfID + "\t" + negWgd + "\t" + df2.format(origLRTH0Full[negWgd-1]) + "\t" + df2.format(medianH0FullDist[negWgd-1]) + "\t" + df2.format(pct95H0FullDist[negWgd-1]) + "\t" + df2.format(FDRpH0Full[k][negWgd-1]) + "\t" + df2.format(origLRTH0Neg[negWgd-1]) + "\t" + df2.format(medianH0NegDist[negWgd-1]) + "\t" + df2.format(pct95H0NegDist[negWgd-1]) + "\t" + df2.format(FDRpH0Neg[k][negWgd-1]) + "\tundecided");
							nrUndecided[k]++;
							nrTopGFsUndecided[negWgd-1]++;
						}
					}
				}
			}

			// Output numbers of WGDs detected/rejected per gene family

			output_perGF.write("GF rank (starts from 0)\tGF ID\t# negative WGDs accepted\t# negative WGDs rejected\t# undecided\t# rejected both H0\n" );

			for (int k=0 ; k<GFs.length; k++) {
				// Start at 0 in arrays and such
				int i = GFs[k]; // i = GF number
				// !! take into account that numbering of GFs starts from 0 for final output, WGD numbering starts at 1
				int gfNumber = i;
				// Get GF ID from output file "i" from Negative observedLR directory
				String fileName = String.format("output_%d.txt",i);
				String origFile = new File(origLRTDir,fileName).toString();
				List<List<String>> outputOrigLRT = cmf.readMapFile(origFile);
				// GF ID is the first field on the first line
				gfID = outputOrigLRT.get(0).get(0);

				output_perGF.write(gfNumber + "\t" + gfID + "\t" + nrNegWGDsAccepted[k] + "\t" + nrNegWGDsRejected[k] + "\t" + nrUndecided[k] + "\t" + nrBothRejected[k] + "\n" );
			}

			// Output how many GFs detect/reject particular WGD

			output_perWGD_Top.write("negative WGD number (starts from 1)\t# Top GFs supporting negative WGD\t# Top GFs rejecting negative WGD\t# Top GFs undecided\t# Top GFs that rejected both H0\n" );

			for (int negWgd = 1; negWgd<=wgdnumber; negWgd++) {
				output_perWGD_Top.write(negWgd + "\t" + nrTopGFsAccept[negWgd-1] + "\t" + nrTopGFsReject[negWgd-1] + "\t" + nrTopGFsUndecided[negWgd-1] + "\t" + nrTopGFsRejectBoth[negWgd-1] + "\n" );
			}

			output_perGF.close();
			output_perWGD_Top.close();
		}
		catch (IOException e) {
			System.out.println(e);
		}
	}
}