package workflows.positiveWGMs;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import utils.bdmodel.EstimateLambdaRootSize;
import utils.bdmodel.GenerateObservations;
import utils.bdmodel.Node;
import utils.bdmodel.ProbCalculator;
import utils.parsers.SpeciesTreeParser;
import utils.parsers.WGMparser;
import utils.parsers.CommonFunctions;

// Class Declaration
public class PositiveWGM_H0Rm {

	// Fields (instance variables): attributes or properties of the class
	private final String treeFile;
	private final String wgdFile;
	private final int gfNumber;
	private final String combinedOutputFile;
	private final String lrtOriginalFile;

	public final static int defmaxNodeSize = 100; 
	public final static int defMinNodeSize = 1;
	public final static double partitionSize = 0.1;
	public final static int lenMCMC = 1000;
	public int maxRootSize = 10;
	public EstimateLambdaRootSize rvsAcc = new EstimateLambdaRootSize(maxRootSize);

	// Constructor to initialize an object of this class
	public PositiveWGM_H0Rm(String treeFile, String wgdFile, int gfNumber, String combinedOutputFile, String lrtOriginalFile) {
		this.treeFile = treeFile;
		this.wgdFile = wgdFile;
		this.gfNumber = gfNumber;
		this.combinedOutputFile = combinedOutputFile;
		this.lrtOriginalFile = lrtOriginalFile;
	}

	// Methods: behaviors or operations that objects of the class can perform
	public static double calculateLRT(double loglkH0, double loglkH1) {
		return (2 * (loglkH1 - loglkH0));
	}

	public double[] calcOptRLamLoglkGF(Node rootOfPartitionedTree, int[] obs, int maxNodeSize, ProbCalculator prc) {
		rootOfPartitionedTree.setLeafValues(obs);
		ArrayList<Node> speTree = SpeciesTreeParser.setMaxNodeSize(rootOfPartitionedTree, maxNodeSize);
		return rvsAcc.cpm_optimize_r_lambda(speTree, rootOfPartitionedTree, prc);
	}

	public void execute() {

		/* ******************************************************************** */
		/* H_0: Current PosWGM (true positive) removed
		/* H_1: All PosWGMs present
		/* ******************************************************************** */

		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat df2 = new DecimalFormat("0.00000");
		ProbCalculator probCalc = new ProbCalculator();
		CommonFunctions cmf = new CommonFunctions();
		
		// Note: WGM file can only contain true WGMs, no negatives
		WGMparser wgm = new WGMparser();
		List<List<String>> wgdList = wgm.readInputFile(wgdFile);
		int posWGDcount= wgdList.size();
		
		int startpos = 0;
		int endpos = posWGDcount-1;


		/* Read gene family ID from file listing the ranked GFs */
		/* **************************************************** */

		// Get GF ID by index (ranking starts from 0)
		List<List<String>> combinedOutput = cmf.readMapFile(combinedOutputFile);
		List<String> myGFrecord = combinedOutput.get(gfNumber);
		String gfID = myGFrecord.get(0);


		/* Generate a simulated gene count profile based on H0Rm (continues later) */
		/* *********************************************************************** */

		// For simulated LRTs: repeat 1000 times = array jobs 
		GenerateObservations go = new GenerateObservations(0, defMinNodeSize, defmaxNodeSize, false, probCalc, lenMCMC);


		/* Print gene family ID */
		/* ******************** */

		System.out.println(gfID);		


		/* Get tree with complete set of positive WGMs */
		/* ******************************************* */

		// Note: WGM file can only contain true WGMs, no negatives, otherwise buildInsertWGDsandPartitionTree does not work
		Node rootFullTree = SpeciesTreeParser.buildInsertWGDsandPartitionTree(treeFile, wgdFile, partitionSize,
				defmaxNodeSize);
		

		/* Loop through the "positive" WGMs (we are sure they occurred) */
		/* ************************************************************ */

		for (int wgd = startpos; wgd <= endpos; wgd++) {

			/* Get tree with current positive WGM removed */
			/* ****************************************** */

			// Code works for removal of one of 2 or 3 WGMs on same branch (max 3!)
			Node rootRmWGD = new Node();
			rootRmWGD = SpeciesTreeParser.buildAndPartitionTree_ReverseEng(treeFile, wgdFile,partitionSize, defmaxNodeSize, wgd);
			
			
			/* Read gene family ID, previously optimized lambda and root size under H0Rm */
			/* ************************************************************************* */

			// From the first row of ObservedLRT output file, get only gene family ID
			// Note: column format in file is [ID, root, lambda, loglik]
			List<List<String>> originalOptParamsRmWGD = cmf.readMapFile(lrtOriginalFile);
			
			String gfID2 = originalOptParamsRmWGD.get(0).get(0);
			if (!gfID.equalsIgnoreCase(gfID2)) {
				System.out.println("error");
			}			
			
			// From the row related to the current WGM, get root size and lambda under the "Rm" scenario
			// Note: column format in file is [root, lambda, loglik, loglik_ratio]
			// Note: to get the correct row, "+1" compensates for the first row of such output containing parameters for the "Full" scenario
			int lineNumberWGD = wgd + 1;
			List<String> lineWGD = originalOptParamsRmWGD.get(lineNumberWGD);
			int rootSizeStarRmWGD = (int) Math.floor(Double.parseDouble(lineWGD.get(0)));
			double RmTree_lambdaStar = Double.parseDouble(lineWGD.get(1));

		
			/* Generate a simulated gene count profile based on H0Rm */
			/* ***************************************************** */

			int[] obs_H0 = go.generateObservation(rootRmWGD, rootSizeStarRmWGD, RmTree_lambdaStar);


			/* Calculate root size, lambda and log-likelihood under H1Full */
			/* *********************************************************** */

			double[] fullTree_simObs_optParams = this.calcOptRLamLoglkGF(rootFullTree, obs_H0, defmaxNodeSize,
					probCalc);


			/* Calculate root size, lambda and log-likelihood under H0Rm */
			/* ********************************************************* */

			double[] rmWGD_sim_optParams = this.calcOptRLamLoglkGF(rootRmWGD, obs_H0, defmaxNodeSize, probCalc);


			/* Print root size, lambda and log-likelihood under H0Rm */
			/* ***************************************************** */

			System.out.print(rmWGD_sim_optParams[0] + "\t" + df.format(rmWGD_sim_optParams[1]) + "\t"
					+ df2.format(rmWGD_sim_optParams[2]) + "\t");
			
			
			/* Print root size, lambda and log-likelihood under H1Full */
			/* ******************************************************* */
			
			System.out.print(fullTree_simObs_optParams[0] + "\t" + df.format(fullTree_simObs_optParams[1]) + "\t"
					+ df2.format(fullTree_simObs_optParams[2]) + "\t");

	
			/* Calculate and print the log-likelihood ratio between H0Rm and H1Full */
			/* ******************************************************************** */

			double lrtSim = calculateLRT(rmWGD_sim_optParams[2], fullTree_simObs_optParams[2]);

			System.out.print(df2.format(lrtSim) + "\n");
		}
	}
}