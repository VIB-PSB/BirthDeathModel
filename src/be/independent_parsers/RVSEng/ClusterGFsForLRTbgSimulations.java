package be.ugent.psb.setas.independent_parsers.RVSEng;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class ClusterGFsForLRTbgSimulations {

	public static void main(String [] args){
		
		String lrtBGfile = "/home/setas/Desktop/setas/Project1/Simulations/LRTbackground/sorted_TauMon2Musa3_RmWGDs.txt";
		CommonFunctions cmf = new CommonFunctions();
		
		ArrayList<Integer> rootSizes= cmf.readColX_Int(lrtBGfile, 1);
		
		List<List<String>> map = cmf.readMapFile(lrtBGfile);
		
		
		for(int i=0; i<9178; i++){
			
			if(rootSizes.get(i) > 5){
				
				List<String> ls = map.get(i);
			
				for(int j=0; j<ls.size();j++){
					
					System.out.print(ls.get(j)+"\t");
				}
				
				System.out.print("\n");
			}

		}
		
		
	}
	
	
}
