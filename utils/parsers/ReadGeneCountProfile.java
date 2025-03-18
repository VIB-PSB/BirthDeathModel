package utils.parsers;

import java.util.List;

public class ReadGeneCountProfile {

	public int[] findObservationBasedOnGFid(String gfcountfile, String gfID, int numberOfLeaves) {

		int[] observation = new int[numberOfLeaves];
		CommonFunctions comFunc = new CommonFunctions();
		List<List<String>> map = comFunc.readMapFile(gfcountfile);

		for (int i = 0; i < map.size(); i++) {
			if (map.get(i).get(0).equals(gfID)) {
				for (int j = 0; j < observation.length; j++) {
					observation[j] = Integer.parseInt(map.get(i).get(j + 3));
				}
				return observation;
			}
		}

		return observation;
	}
	
}
