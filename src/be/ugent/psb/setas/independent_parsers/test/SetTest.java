package be.ugent.psb.setas.independent_parsers.test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class SetTest {

	private void doTest() {
		HashSet<String> a = new HashSet<String>();
		HashSet<String> b = new HashSet<String>();

		a.add("A");
		a.add("B");

		b.add("B");
		b.add("A");

		// System.out.println(a.equals(b));

		HashMap<Set<String>, String> testMap = new HashMap<Set<String>, String>();
		testMap.put(a, "Tadaaaa");

		// System.out.println(testMap.get(b));
	}

	public static void main(String[] args) {
		// SetTest test = new SetTest();
		//
		// test.doTest();

		CommonFunctions cmf = new CommonFunctions();
		String myPath = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/MusaFirst2WGDs/basedOnLambdaRanking/TopBottom/TopBottom_lamRank_MusaFirst2Close_ExpFunctPlots";

		List<List<String>> map = cmf.readMapFile(myPath);

		for (List<String> lsm : map) {
			
			if (!lsm.get(1).equalsIgnoreCase("ORTHO0002576_1") && !lsm.get(1).equalsIgnoreCase("ORTHO002324")
					&& !lsm.get(1).equalsIgnoreCase("ORTHO010110")) {

				for (int i = 0; i < lsm.size() - 1; i++) {

					System.out.print(lsm.get(i) + "\t");
				}

				System.out.print(lsm.get(lsm.size() - 1) + "\n");
			}
		}

	}
}
