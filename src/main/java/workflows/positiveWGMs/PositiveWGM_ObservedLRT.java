package workflows.positiveWGMs;

import utils.parsers.ReadGeneCountProfile;
import java.util.ArrayList;
import java.util.List;
import utils.bdmodel.Node;
import utils.bdmodel.ProbCalculator;
import utils.bdmodel.EstimateLambdaRootSize;
import utils.parsers.SpeciesTreeParser;
import utils.parsers.WGMparser;
import utils.parsers.CommonFunctions;

// Class Declaration
public class PositiveWGM_ObservedLRT {

	// Fields (instance variables): attributes or properties of the class
	private final String treeFile;
	private final String wgdFile;
	private final int gfNumber;
        private final int testedWGD;
	private final String combinedOutputFile;
	private final String gfCountsFile;

	public final static int defmaxNodeSize = 100;
	public final static int defMinNodeSize = 1;
	public final static double partitionSize = 0.1;
	public final static int lenMCMC = 5000;
	public int maxRootSize = 10;
	public EstimateLambdaRootSize rvsAcc = new EstimateLambdaRootSize(maxRootSize);

	// Constructor to initialize an object of this class
	public PositiveWGM_ObservedLRT(String treeFile, String wgdFile, int testedWGD, int gfNumber, String combinedOutputFile, String gfCountsFile) {
		this.treeFile = treeFile;
		this.wgdFile = wgdFile;
                this.testedWGD = testedWGD;
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

		/* ****************************************************************** */
		/* H_0: All PosWGMs (true positive) present
		/* H_1: Current PosWGM removed
		/* ****************************************************************** */

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

		List<List<String>> combinedOutput = cmf.readMapFile(combinedOutputFile);
		List<String> myGFrecord = combinedOutput.get(gfNumber);
		String gfID = myGFrecord.get(0);
		// int rootSizeStar = Integer.parseInt(myGFrecord.get(1)); // Not used
		// double fullTree_lambdaStar = Double.parseDouble(myGFrecord.get(2)); // Not used


		/* Get tree with complete set of positive WGMs */
		/* ******************************************* */

		Node rootFullTree = SpeciesTreeParser.buildInsertWGDsandPartitionTree(treeFile, wgdFile, partitionSize,
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
		

		/* Loop through the "positive" WGMs (we are sure they occurred) */
		/* ************************************************************ */

                if(testedWGD == -1){
                    
                    for (int posWgd = startpos; posWgd <= endpos; posWgd++) {
			
			/* Get tree with current positive WGM removed */
			/* ****************************************** */

			// Code supports removal of one of 2 or 3 WGMs on same branch (max 3!)
			Node rootRmWGD = new Node();	
			rootRmWGD = SpeciesTreeParser.buildAndPartitionTree_ReverseEng(treeFile, wgdFile,partitionSize, defmaxNodeSize, posWgd);


			/* Calculate and print root size, lambda and log-likelihood under H1Rm */
			/* ******************************************************************* */

			double[] posWGD_orgObs_optParams = this.calcOptRLamLoglkGF(rootRmWGD, originalGFCounts,
					defmaxNodeSize, probCalc);

			System.out.print(posWGD_orgObs_optParams[0] + "\t" + posWGD_orgObs_optParams[1] + "\t"
					+ posWGD_orgObs_optParams[2] + "\t");
			

			/* Calculate and print the log-likelihood ratio between H0Full and H1Rm */
			/* ******************************************************************** */

			double lrtOriginal = calculateLRT(fullTree_loglkStar[2], posWGD_orgObs_optParams[2]);			
			
			System.out.print(lrtOriginal + "\n");
                    }
		}
                else{
                    for (int posWgd = startpos; posWgd <= endpos; posWgd++) {
                        
                        if(posWgd==testedWGD){
			
                            /* Get tree with current positive WGM removed */
                            /* ****************************************** */

                            // Code supports removal of one of 2 or 3 WGMs on same branch (max 3!)
                            Node rootRmWGD = new Node();	
                            rootRmWGD = SpeciesTreeParser.buildAndPartitionTree_ReverseEng(treeFile, wgdFile,partitionSize, defmaxNodeSize, posWgd);


                            /* Calculate and print root size, lambda and log-likelihood under H1Rm */
                            /* ******************************************************************* */

                            double[] posWGD_orgObs_optParams = this.calcOptRLamLoglkGF(rootRmWGD, originalGFCounts,
                                            defmaxNodeSize, probCalc);

                            System.out.print(posWGD_orgObs_optParams[0] + "\t" + posWGD_orgObs_optParams[1] + "\t"
                                            + posWGD_orgObs_optParams[2] + "\t");


                            /* Calculate and print the log-likelihood ratio between H0Full and H1Rm */
                            /* ******************************************************************** */

                            double lrtOriginal = calculateLRT(fullTree_loglkStar[2], posWGD_orgObs_optParams[2]);			

                            System.out.print(lrtOriginal + "\n");
                        }
                        else{
                            System.out.print("0\t0\t0\t0\n");
                        }
                    }

                }
	}
}