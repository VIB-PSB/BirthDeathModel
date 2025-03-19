package be.ugent.psb.setas.independent_parsers.GO;

import java.util.ArrayList;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class ParsePlazaGoAnnotationFile_BiNGO {

	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();

		String GFidGeneIDfile = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/AllGenes_AllGFs_37pePlazaIDs";

		ArrayList<String> myGFID = cmf.readColX_String(GFidGeneIDfile, 0);
		ArrayList<String> myGeneID = cmf.readColX_String(GFidGeneIDfile, 1);
		
		String file2="/home/setas/Desktop/AnnotationPlaza/go.zma.csv";
		ArrayList<String> PlazaGeneID = cmf.readColX_String_SemiColon(file2, 2);
		ArrayList<String> PlazaGOID= cmf.readColX_String_SemiColon(file2,3);
		
//		System.out.println(PlazaGOID.size());
		
		ArrayList<String> PlazaGOID_refined = new ArrayList<String>(PlazaGOID.size());
		
		for(int i=0; i<PlazaGOID.size();i++){
			
//			System.out.println(PlazaGOID.get(i).split(":")[1]);
			
			PlazaGOID_refined.add(i, PlazaGOID.get(i).split(":")[1]);
		}        
        
		
		String tag= "ZM";
        
        for(int j=0; j<myGFID.size();j++){
        	
        	String currentGFid = myGFID.get(j);
        	String currentGeneID = myGeneID.get(j);
        
        	//&& currentGeneID.split("ATR").length <= 1 //only for Athaliana
        	
        	if(currentGeneID.split(tag).length >1 ){
        	for(int i=0; i<PlazaGeneID.size();i++){
        		
//        		System.out.println(currentGeneID+"\t"+PlazaGeneID.get(i));
        		
        		if(PlazaGeneID.get(i).equals(currentGeneID)){
        			
        			System.out.println(currentGFid+"="+PlazaGOID_refined.get(i));
        		}
        		
        	}
        	}
        	
        }
        
		
	}

}
