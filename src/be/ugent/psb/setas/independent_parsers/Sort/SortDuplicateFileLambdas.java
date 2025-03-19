package be.ugent.psb.setas.independent_parsers.Sort;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class SortDuplicateFileLambdas {

	public static void main(String[] args) {

		CommonFunctions cmmFunct = new CommonFunctions();

		String path1 = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/9178coreGF_inOrderLam";
		List<String> GF_id_order = cmmFunct.read1ColFile_String(path1);

		String path2 = "/home/setas/git/IndependentParsers/Independent Parsers/src/files/genefamily.angiosperm.dup.rootsupport70.paml.Atha.txt";
		ArrayList<String> GF_id = cmmFunct.readColX_String(path2, 0);
		List<String> Dup_node = cmmFunct.readColX_String(path2, 1);
		List<String> gene_1 = cmmFunct.readColX_String(path2, 2);
		List<String> gene_2 = cmmFunct.readColX_String(path2, 3);
		List<Double> t = cmmFunct.readColX_double(path2, 4);
		List<Double> S = cmmFunct.readColX_double(path2, 5);
		List<Double> N = cmmFunct.readColX_double(path2, 6);
		List<Double> omega = cmmFunct.readColX_double(path2, 7);

		List<Double> dN = cmmFunct.readColX_double(path2, 8);
		List<Double> dS = cmmFunct.readColX_double(path2, 9);

		// for (String gfId_prob : GF_id_order) {
		for (int j = 0; j < 709; j++) {

			String gfId_prob = GF_id_order.get(j);

			if (!cmmFunct.searchListString_boolean(gfId_prob, GF_id)) {
				// System.out.println(gfId_prob+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+"NA");
			}

			else {
				for (int i = 0; i < GF_id.size(); i++) {

					if (GF_id.get(i).equals(gfId_prob)) {

						System.out.print(gfId_prob + "\t" + Dup_node.get(i)
								+ "\t" + gene_1.get(i).split("\\|")[1] + "\t"
								+ gene_2.get(i).split("\\|")[1] + "\t"
								+ t.get(i) + "\t" + S.get(i) + "\t" + N.get(i)
								+ "\t" + omega.get(i) + "\t" + dN.get(i) + "\t"
								+ dS.get(i));
						System.out.println();
					}
				}

			}

		}

	}
}