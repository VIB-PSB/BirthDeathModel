package be.ugent.psb.setas.independent_parsers.KnKs;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class CountNumberOfpairsInWGDRange {

	public static void main(String[] args) {

		CommonFunctions cmm = new CommonFunctions();

		String path = "/mnt/shares/biocomp/groups/group_esb/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/KnKs/GFsOrder_genDup/Ks5_minKs/TopBottom/GenPairs_PAML/Atha_Ks5_Top709.txt";

		double maxrange = 1.1;
		double minrange = 0.5;

		ArrayList<String> gfId = cmm.readColX_String(path, 0);
		ArrayList<String> gene1 = cmm.readColX_String(path, 1);
		ArrayList<String> gene2 = cmm.readColX_String(path, 2);
		ArrayList<Double> t = cmm.readColX_double(path, 3);
		ArrayList<Double> s = cmm.readColX_double(path, 4);
		ArrayList<Double> n = cmm.readColX_double(path, 5);
		ArrayList<Double> omega = cmm.readColX_double(path, 6);
		ArrayList<Double> Kn = cmm.readColX_double(path, 7);
		ArrayList<Double> Ks = cmm.readColX_double(path, 8);

		int numOfGFFromTheTop = 500;
		int counter = 0;
		int i = 0;
		int numOfPairsWithKsInrange = 0;
		
		double sumKnInRange =0;

		String gfIDcurrent = gfId.get(0);
		if (Ks.get(0) >= minrange && Ks.get(0) <= maxrange) {
			numOfPairsWithKsInrange += 1;
		}

		while (counter < numOfGFFromTheTop) {

			i++;
			String gf = gfId.get(i);

			if (Ks.get(i) >= minrange && Ks.get(i) <= maxrange) {

				numOfPairsWithKsInrange += 1;
				
				sumKnInRange += Kn.get(i);
			}

			if (!gf.equals(gfIDcurrent)) {

				counter++;

			}

		}
		
		double averageKnInRange = sumKnInRange/numOfPairsWithKsInrange;
		
		System.out.println(numOfPairsWithKsInrange+"\t"+averageKnInRange);

	}

}
