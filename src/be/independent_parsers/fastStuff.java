package be.ugent.psb.setas.independent_parsers;

import java.util.List;

public class fastStuff {

	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();
		String path = "/home/setas/Desktop/GO_enrichments_review/rank_GF_GO_Desc_MusaFirst2Close_CombinedRanking_correct_sorted.txt";

		List<List<String>> map = cmf.readMapFile(path);

		List<String> firstLine = map.get(0);
		String gfID = firstLine.get(1);
		int lineNumber = 1;

		System.out.print(lineNumber + "\t");

		for (int i = 0; i < firstLine.size() - 1; i++) {
			System.out.print(firstLine.get(i) + "\t");
		}
		System.out.print(firstLine.get(firstLine.size() - 1) + "\n");

		for (int i = 1; i < map.size(); i++) {

			List<String> record = map.get(i);
			String gfID_current = record.get(1);

			if (gfID_current.equalsIgnoreCase(gfID)) {

				System.out.print(lineNumber + "\t");

				for (int j = 0; j < record.size() - 1; j++) {
					System.out.print(record.get(j) + "\t");
				}

				System.out.print(record.get(record.size() - 1) + "\n");
			}

			else {

				gfID = gfID_current;
				lineNumber = lineNumber + 1;

				System.out.print(lineNumber + "\t");

				for (int j = 0; j < record.size() - 1; j++) {
					System.out.print(record.get(j) + "\t");
				}

				System.out.print(record.get(record.size() - 1) + "\n");

			}

		}

	}

}
