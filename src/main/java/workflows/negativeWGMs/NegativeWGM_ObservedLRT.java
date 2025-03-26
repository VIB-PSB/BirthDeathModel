package workflows.negativeWGMs;

import utils.parsers.ReadGeneCountProfile;
import java.util.ArrayList;
import java.util.List;

import utils.bdmodel.EstimateLambdaRootSize;
import utils.bdmodel.Node;
import utils.bdmodel.ProbCalculator;
import utils.parsers.SpeciesTreeParser;
import utils.parsers.WGMparser;
import utils.parsers.CommonFunctions;

// Class Declaration
public class NegativeWGM_ObservedLRT {

	// Fields (instance variables): attributes or properties of the class
	private final String treeFile;
	private final String wgdFile;
	private final int gfNumber;
	private final String combinedOutputFile;
	private final String gfCountsFile;

	public final static int defmaxNodeSize = 100;
	public final static int defMinNodeSize = 1;
	public final static double partitionSize = 0.1;
	public final static int lenMCMC = 5000;
	public int maxRootSize = 10;
	public EstimateLambdaRootSize rvsAcc = new EstimateLambdaRootSize(maxRootSize);

	// Constructor to initialize an object of this class
	public NegativeWGM_ObservedLRT(String treeFile, String wgdFile, int gfNumber, String combinedOutputFile, String gfCountsFile) {
		this.treeFile = treeFile;
		this.wgdFile = wgdFile;
		this.gfNumber = gfNumber;
		this.combinedOutputFile = combinedOutputFile;
		this.gfCountsFile = gfCountsFile;
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

		/* **************************************************************************************** */
		/* H_0: Only all PosWGMs (true positive) present
		/* H_1: Current NegWGM added (true negative, i.e. additional WGM which is not really there)
		/* **************************************************************************************** */

		ProbCalculator probCalc = new ProbCalculator();
		CommonFunctions cmf = new CommonFunctions();

		WGMparser wgm = new WGMparser();
		List<List<String>> wgdList = wgm.readInputFile(wgdFile);
		int posWGDcount=0;
		for (int i = 0; i < wgdList.size(); i++) {
			if (wgdList.get(i).get(0).equalsIgnoreCase("negatives")) {
				posWGDcount=i;
			}
		}
		

		// Skip negatives keyword
		int startneg = posWGDcount+1;
		int endneg = wgdList.size()-1;


		/* Read gene family ID from file listing the ranked GFs */
		/* **************************************************** */

		// Get GF ID by index (ranking starts from 0)
		// Note: column format in file is [ID, root, lambda, loglik]
		List<List<String>> combinedOutput = cmf.readMapFile(combinedOutputFile);
		List<String> myGFrecord = combinedOutput.get(gfNumber);
		String gfID = myGFrecord.get(0);
		// int rootSizeStar = Integer.parseInt(myGFrecord.get(1)); // Not used
		// double fullTree_lambdaStar = Double.parseDouble(myGFrecord.get(2)); // Not used


		/* Get tree with complete set of positive WGMs */
		/* ******************************************* */

		Node rootFullTree = SpeciesTreeParser.buildInsertWGDsandPartitionTree_Negatives(treeFile, wgdFile, partitionSize,
				defmaxNodeSize);
		int nrspecies = rootFullTree.getLeaves().size();


		/* Read observed gene count profile */
		/* ******************************** */

		ReadGeneCountProfile calcLrt = new ReadGeneCountProfile();
		int[] originalGFCounts = calcLrt.findObservationBasedOnGFid(gfCountsFile, gfID, nrspecies);
		

		/* Calculate root size, lambda and log-likelihood under H0Full */
		/* *********************************************************** */

		double [] fullTree_loglkStar = this.calcOptRLamLoglkGF(rootFullTree, originalGFCounts,
				defmaxNodeSize, probCalc);

		
		/* Print gene family ID, root size, lambda and log-likelihood under H0Full */
		/* *********************************************************************** */

		System.out.println(gfID + "\t" + fullTree_loglkStar[0] + "\t" + fullTree_loglkStar[1] + "\t" + fullTree_loglkStar[2]);		


		/* Loop through the "negative" WGMs (we are sure they DIDN'T occur) */
		/* **************************************************************** */

		for (int negWgd = startneg; negWgd <= endneg; negWgd++) {	
			
			/* Get tree with current negative WGM added */
			/* **************************************** */

			Node rootNegWGD = new Node();	
			rootNegWGD = SpeciesTreeParser.buildAndPartitionTree_ReverseEng_Negatives(treeFile, wgdFile, partitionSize, defmaxNodeSize, negWgd);
			

			/* Calculate and print root size, lambda and log-likelihood under H1Neg */
			/* ******************************************************************** */

			double[] negWGD_orgObs_optParams = this.calcOptRLamLoglkGF(rootNegWGD, originalGFCounts,
					defmaxNodeSize, probCalc);
			
			System.out.print(negWGD_orgObs_optParams[0] + "\t" + negWGD_orgObs_optParams[1] + "\t"
					+ negWGD_orgObs_optParams[2] + "\t");
	
		
			/* Calculate and print the log-likelihood ratio between H0Full and H1Neg */
			/* ********************************************************************* */

			double lrtOriginal = calculateLRT(fullTree_loglkStar[2], negWGD_orgObs_optParams[2]);

			System.out.print(lrtOriginal + "\n");
		}
	}
}