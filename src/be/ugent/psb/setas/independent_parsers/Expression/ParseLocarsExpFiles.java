package be.ugent.psb.setas.independent_parsers.Expression;

import java.util.ArrayList;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class ParseLocarsExpFiles {
	
	public static void main(String [] args){
				
		CommonFunctions cmf = new CommonFunctions();
		
		String cornetFile ="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/newfromLor/CornetBYadj.txt";
		ArrayList<String> gene1_cornet = cmf.readColX_String(cornetFile, 0);
		ArrayList<String> gene2_cornet = cmf.readColX_String(cornetFile, 1);
		ArrayList<Double> corrolation = cmf.readColX_double(cornetFile, 2);	
		
		String codemlFile ="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/newfromLor/genefamily.angiosperm.core.txt.codeml";
		
		ArrayList<String> gfIDs = cmf.readColX_String(codemlFile, 0);
		ArrayList<String> geneID1 = cmf.readColX_String(codemlFile, 2);
		ArrayList<String> geneID2 = cmf.readColX_String(codemlFile, 3);
		
//		ArrayList<Double> dn = cmf.readColX_double(codemlFile,8 );
		ArrayList<Double> ds = cmf.readColX_double(codemlFile,9);		
		
		for(int i=0; i<gfIDs.size(); i++){			
			String gfIDcurrent = gfIDs.get(i);			
			String gene1 = geneID1.get(i).split("|")[1];
			String gene2 = geneID1.get(i).split("|")[1];
					
			for(String g: gene1_cornet){}
		}			
	}
}
