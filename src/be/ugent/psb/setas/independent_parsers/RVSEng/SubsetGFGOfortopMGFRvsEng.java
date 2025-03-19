package be.ugent.psb.setas.independent_parsers.RVSEng;

import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class SubsetGFGOfortopMGFRvsEng {

	
	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();
	    String path="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/Tair/GF_GOTair_forBingo_orderlam";
		
	    List<String> GFids = cmf.readColX_String_Delimiter(path, 0,"=",true);
	    List<Integer> GOids = cmf.readColX_Int_Delimiter(path, 1,"=",true);
	    
	    String path2="/home/setas/Desktop/setas/Project1/RvsEng/TopGFs/AngioSperms/GOenrichment/top709MarkerGFid_LowerAvgLogLkRm1WGD.txt";
	    List<String> newGFids = cmf.read1ColFile_String(path2);
	    
	
	    
	    for(int i=0 ; i<newGFids.size(); i++){
	    	
	    	String gfProb = newGFids.get(i);
	    	
	    	
	    	for(int j=0 ; j<GFids.size();j++){
	    		
	    		if(GFids.get(j).equals(gfProb)){
	    			
	    			System.out.println(gfProb+"="+GOids.get(j));
	    		}
	    	}
	    	
	    }
	
	
	}
	
}
