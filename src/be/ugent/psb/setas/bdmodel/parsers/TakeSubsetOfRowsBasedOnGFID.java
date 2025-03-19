package be.ugent.psb.setas.bdmodel.parsers;

import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class TakeSubsetOfRowsBasedOnGFID {

	
	public static void main(String args[]) {
		
		CommonFunctions cmf = new CommonFunctions();
		String path_combOutFile = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/MusaFirst2WGDs/CombinedOutput/combinedOutput_MusaFirst2close";
		String path_selectedGFIDs = "/home/setas/Desktop/GFs_selectedJackKnife";
		
		List<List<String>> map = cmf.readMapFile(path_combOutFile);
		List<String> gfIDs = cmf.read1ColFile_String(path_selectedGFIDs);
		
		for(int i=0; i<map.size();i++) {
			
			String gfID_comb = map.get(i).get(0);
			List<String> line = map.get(i);
			int lenxLine = line.size();
			
			if(gfIDs.contains(gfID_comb)) {
				
				for(int j=0; j< lenxLine-1 ;j++) {
					
					System.out.print(line.get(j)+"\t");
				}
				
				System.out.print(line.get(lenxLine-1)+"\n");
			}
		}
		
	}
}
