package be.ugent.psb.setas.bdmodel.model.RVS_Engineering;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import be.ugent.psb.setas.bdmodel.model.GenerateObservations;
import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.bdmodel.model.ProbCalculator;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;
import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class LRTbackgroundDist_RvsEng {
	/**
	 * method 1: read the mgf file, pick mgf number i, set the root size and //
	 * lambda according to that, // generate one observation based on that, set the
	 * leaf sizes, calculate the // loglks under original and rmWG0 model, //
	 * calculate the Q(observation) = 2(Loglk_H0 - Loglk_H1) // repeat it for 1000
	 * observations // count how many times Q-original > Q(observation-i) = quantile
	 * // if quantile >= 95% , reject H0
	 **/
	public ProbCalculator probCalc;
	public final static int defmaxNodeSize = 100; // min should be 23.
	public final static int defMinNodeSize = 1;
	public final static double partitionSize = 0.1;
	public final static int lenMCMC = 5000;
	public int maxRootSize = 10;
	public RVSEng_accurateLambda_MusaFirst2 rvsAcc = new RVSEng_accurateLambda_MusaFirst2(maxRootSize);

	public LRTbackgroundDist_RvsEng(ProbCalculator probCalc) {
		this.probCalc = probCalc;
	}

	public static double calculateLRT(double loglkH0, double loglkH1) {
		return (2 * (loglkH1 - loglkH0));
	}

	public double[] calcOptRLamLoglkGF(Node rootOfPartitionedTree, int[] obs, int maxNodeSize, ProbCalculator prc) {

		rootOfPartitionedTree.setLeafValues(obs);// It does set the values of leaves correctly, in postOrder manner, as
													// the obs array was generated.
		ArrayList<Node> speTree = SpeciesTreeParser.setMaxNodeSize(rootOfPartitionedTree, maxNodeSize);
		return rvsAcc.cpm_optimize_r_lambda(speTree, rootOfPartitionedTree, prc);
	}

	public static void main(String[] args) {

		/** H0: all WGDs (A), H1: RmWGDi (B) */
		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat df2 = new DecimalFormat("0.00000");
		ProbCalculator probCalc = new ProbCalculator();
		CommonFunctions cmf = new CommonFunctions();
		String treeFile = args[0];
		String wgdFile = args[1];
		int gfNumber = Integer.parseInt(args[2]); // manually set to 0 for now, later I submit 200 such scripts?!
		String wgdFile_rmWGD8 = args[3];
		String wgdFile_rmWGD9 = args[4];
		String wgdFile_rmWGD12 = args[5];
		String wgdFile_rmWGD16 = args[6];
		String combinedOutputFile = args[7];

//		String treeFile = "C:\\Users\\setar\\git\\StochasticBD\\37spe-MGCF5.txt";
//		String wgdFile = "C:\\Users\\setar\\git\\StochasticBD\\WGDs-MusaFirst2Successive";
//		int gfNumber = 5;
//		String wgdFile_rmWGD8 = "C:\\Users\\setar\\git\\StochasticBD\\WGDs-MusaFirst2Successive_rm8";
//		String wgdFile_rmWGD9 = "C:\\Users\\setar\\git\\StochasticBD\\WGDs-MusaFirst2Successive_rm9";
//		String wgdFile_rmWGD12 = "C:\\Users\\setar\\git\\StochasticBD\\WGDs-MusaFirst2Successive_rm12";
//		String wgdFile_rmWGD16 = "C:\\Users\\setar\\git\\StochasticBD\\WGDs-MusaFirst2Successive_rm16";
//		String combinedOutputFile = "C:\\Users\\setar\\git\\StochasticBD\\combinedOutput_MusaFirst2close_orderedLambda";

		LRTbackgroundDist_RvsEng lrtDist = new LRTbackgroundDist_RvsEng(probCalc);

		List<List<String>> combinedOutput = cmf.readMapFile(combinedOutputFile);
		List<String> myGFrecord = combinedOutput.get(gfNumber);
		String gfID = myGFrecord.get(0);
//		int rootSizeStar = Integer.parseInt(myGFrecord.get(1));
//		double fullTree_lambdaStar = Double.parseDouble(myGFrecord.get(2));

		/* For original LRTs */
//		CalculateCorrectLRToriginalFiles calcLrt = new CalculateCorrectLRToriginalFiles();
//		//String gfCountsFile = "C:\\Users\\setar\\git\\StochasticBD\\37spe-MGCF5-9178core-OrderNewickTree-noHeader";
//		String gfCountsFile = args[8];
//		
//		Node rootFullTree = SpeciesTreeParser.buildInsertWGDsandPartitionTree(treeFile, wgdFile, partitionSize,
//				defmaxNodeSize);
//		int[] originalGFCounts = calcLrt.findObservationBasedOnGFid(gfCountsFile, gfID, 37);
//		
//		double [] fullTree_loglkStar = lrtDist.calcOptRLamLoglkGF(rootFullTree, originalGFCounts,
//				defmaxNodeSize, probCalc);		
//		System.out.println(gfID + "\t" + fullTree_loglkStar[0] + "\t" + fullTree_loglkStar[1] + "\t" + fullTree_loglkStar[2]);		
//		
//		for (int wgd = 0; wgd < 20; wgd++) {			
//			Node rootRmWGD = new Node();		
//			rootRmWGD = SpeciesTreeParser.buildAndPartitionTree_RVS(treeFile, wgdFile, wgdFile_rmWGD8, wgdFile_rmWGD9,
//					wgdFile_rmWGD12, wgdFile_rmWGD16, partitionSize, defmaxNodeSize, wgd);
//			
//			double[] rmWGD_orgObs_optParams = lrtDist.calcOptRLamLoglkGF(rootRmWGD, originalGFCounts,
//					defmaxNodeSize, probCalc);		
//			System.out.print(rmWGD_orgObs_optParams[0] + "\t" + rmWGD_orgObs_optParams[1] + "\t"
//					+ rmWGD_orgObs_optParams[2] + "\t");
//			
//			double lrtOriginal = calculateLRT(fullTree_loglkStar[2], rmWGD_orgObs_optParams[2]);			
//			System.out.print(lrtOriginal + "\n");}

//		/** For simulated LRTs: repeat 1000 times = array jobs H0 = FullTree*/
//		GenerateObservations go = new GenerateObservations(0, defMinNodeSize, defmaxNodeSize, false, probCalc, lenMCMC);
//
//		Node rootFullTree_noPartition = SpeciesTreeParser.buildInsertWGDsandPartitionTree(treeFile, wgdFile, 0,
//				defmaxNodeSize); // build without partition for generating obs
//
//		int[] obs_H0 = go.generateObservation(rootFullTree_noPartition, rootSizeStar, fullTree_lambdaStar);
//
//		System.out.println(gfID);
//		
//		for (int o : obs_H0) {System.out.print(o + "\t");}
//		System.out.print("\n");
//
//		// rebuild and insert partitions for keeping lambda in range:
//		Node rootFullTree = SpeciesTreeParser.buildInsertWGDsandPartitionTree(treeFile, wgdFile, partitionSize,
//				defmaxNodeSize);
//
//		double[] fullTree_simObs_optParams = lrtDist.calcOptRLamLoglkGF(rootFullTree, obs_H0, defmaxNodeSize, probCalc);
//
//		System.out.println(fullTree_simObs_optParams[0] + "\t" + df.format(fullTree_simObs_optParams[1]) + "\t"
//				+ df2.format(fullTree_simObs_optParams[2]));
//
//		for (int wgd = 0; wgd < 20; wgd++) {
//
//			Node rootRmWGD = new Node();
//			rootRmWGD = SpeciesTreeParser.buildAndPartitionTree_RVS(treeFile, wgdFile, wgdFile_rmWGD8, wgdFile_rmWGD9,
//					wgdFile_rmWGD12, wgdFile_rmWGD16, partitionSize, defmaxNodeSize, wgd);
//
//			double[] rmWGD_sim_optParams = lrtDist.calcOptRLamLoglkGF(rootRmWGD, obs_H0, defmaxNodeSize, probCalc);
//
//			System.out.print(rmWGD_sim_optParams[0] + "\t" + df.format(rmWGD_sim_optParams[1]) + "\t"
//					+ df2.format(rmWGD_sim_optParams[2]) + "\t");
//			double lrtSim = calculateLRT(fullTree_simObs_optParams[2], rmWGD_sim_optParams[2]);
//
//			System.out.print(df2.format(lrtSim) + "\n");
//		}

		/** For simulated LRTs: repeat 1000 times = array jobs H0=Rm-WGD */
		//CalculateCorrectLRToriginalFiles calcLrt = new CalculateCorrectLRToriginalFiles();
		GenerateObservations go = new GenerateObservations(0, defMinNodeSize, defmaxNodeSize, false, probCalc, lenMCMC);

		System.out.println(gfID);
		// output_0 corresponding to the top gene family is missing.. we have to do it
		// manually later.
//		String pathOrgLrtDir = "C:\\Users\\setar\\git\\StochasticBD\\src\\OriginalLRTs_H0Full_Top100";	
//		String outputFileName = "\\output_"; //For Windows
//		String outputFileName = "/output_"; //For linux servers
//		String gfNumber_string = String.valueOf(gfNumber);
//		String outputfilePath = outputFileName.concat(gfNumber_string);	

//		String lrtOriginalFile_H0Full = "C:\\Users\\setar\\git\\StochasticBD\\src\\OriginalLRTs_H0Full_Top100\\output_5";
		String lrtOriginalFile_H0Full = args[8];

		Node rootFullTree = SpeciesTreeParser.buildInsertWGDsandPartitionTree(treeFile, wgdFile, partitionSize,
				defmaxNodeSize);

		for (int wgd = 0; wgd < 20; wgd++) {

			Node rootRmWGD = new Node();
			rootRmWGD = SpeciesTreeParser.buildAndPartitionTree_RVS(treeFile, wgdFile, wgdFile_rmWGD8, wgdFile_rmWGD9,
					wgdFile_rmWGD12, wgdFile_rmWGD16, partitionSize, defmaxNodeSize, wgd);

			//because first line of such outputs is the GFid and so on..
			int lineNumberWGD = wgd + 1;

			List<List<String>> originalOptParamsRmWGD = cmf.readMapFile(lrtOriginalFile_H0Full);
			
			String gfID2 = originalOptParamsRmWGD.get(0).get(0);
			
			if (!gfID.equalsIgnoreCase(gfID2)) {
				System.out.println("error");
			}
			List<String> lineWGD = originalOptParamsRmWGD.get(lineNumberWGD);
			int rootSizeStarRmWGD = (int) Math.floor(Double.parseDouble(lineWGD.get(0)));
			double RmTree_lambdaStar = Double.parseDouble(lineWGD.get(1));
			//double RmTree_loglkStar = Double.parseDouble(lineWGD.get(2));

			int[] obs_H0 = go.generateObservation(rootRmWGD, rootSizeStarRmWGD, RmTree_lambdaStar);

			double[] fullTree_simObs_optParams = lrtDist.calcOptRLamLoglkGF(rootFullTree, obs_H0, defmaxNodeSize,
					probCalc);
			double[] rmWGD_sim_optParams = lrtDist.calcOptRLamLoglkGF(rootRmWGD, obs_H0, defmaxNodeSize, probCalc);

			System.out.print(rmWGD_sim_optParams[0] + "\t" + df.format(rmWGD_sim_optParams[1]) + "\t"
					+ df2.format(rmWGD_sim_optParams[2]) + "\t");

			System.out.print(fullTree_simObs_optParams[0] + "\t" + df.format(fullTree_simObs_optParams[1]) + "\t"
					+ df2.format(fullTree_simObs_optParams[2]) + "\t");

			double lrtSim = calculateLRT(rmWGD_sim_optParams[2], fullTree_simObs_optParams[2]);

			System.out.print(df2.format(lrtSim) + "\n");
		}
	}
}