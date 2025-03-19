package be.ugent.psb.setas.independent_parsers.TandemBlock;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class PrintTBpercentages_mibelFile {
	
	
	public static void main(String[] args) {

	CommonFunctions cmmFunct = new CommonFunctions();
	
//	String path = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/Lusi/geneIDSpeTB_Lus.txt";
	String path = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/Michiel/geneIDSpeTB.txt";

	ArrayList<String> geneIDs = cmmFunct.readColX_String(path, 0);
	ArrayList<Integer> isBlock = cmmFunct.readColX_Int(path, 2);
	ArrayList<Integer> isTandem = cmmFunct.readColX_Int(path, 3);

//	String path1 = "/home/setas/Desktop/AllGenes_29spePlazaID";
//	ArrayList<String> mygfID = cmmFunct.readColX_String(path1, 0);
//	ArrayList<String> mygeneID = cmmFunct.readColX_String(path1,1);
	
//	String path1 = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/Lusi/AllGenes_29spePlazaID_Lus";
	String path1 = "/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/new/AllGenes_AllGFs_37pePlazaIDs";
	ArrayList<String> mygfID = cmmFunct.readColX_String(path1, 0);
	ArrayList<String> mygeneID = cmmFunct.readColX_String(path1,1);
	
	String [] newTags = new String [7];
	newTags[0]="MA";
	newTags[1]="SI";
	newTags[2]="CM";
	newTags[3]="Csa";
	newTags[4]="HV";
	newTags[5]="GR";
	newTags[6]="TC";
	
	
	for(int i=0; i<mygfID.size();i++){
//	for(int i=0; i<1000;i++){
		
		
		for (String tag : newTags){
		if(mygeneID.get(i).split(tag).length > 1){
		
		int index = cmmFunct.searchListString_index(mygeneID.get(i), geneIDs);
		
		if(index >= 0){
			
			int TB =0;
			if(isBlock.get(index)==1 && isTandem.get(index)==1){
				
				TB=1;
			}
			
			System.out.println(mygfID.get(i)+"\t"+mygeneID.get(i)+"\t"+isBlock.get(index)+"\t"+isTandem.get(index)+"\t"+TB);
		}
	}
	}
		
	}
	
	}
}
