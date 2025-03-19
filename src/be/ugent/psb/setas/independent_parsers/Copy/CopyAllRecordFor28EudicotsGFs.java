package be.ugent.psb.setas.independent_parsers.Copy;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class CopyAllRecordFor28EudicotsGFs {

	
	
	public static void main(String[] args) {
		
		CommonFunctions cmmFun = new CommonFunctions();
		CopyAllRecordFor28EudicotsGFs cpAllRec = new CopyAllRecordFor28EudicotsGFs();
		
//		ArrayList<String> EudGFsInOrder = cmmFun.read1ColFile_String("/home/setas/git/IndependentParsers/Independent Parsers/src/files/GFsOrder_28Eud_37spe_9178coreGF");
////		List<List<String>> map = cmmFun.readMapFile("/home/setas/git/IndependentParsers/Independent Parsers/src/files/GFid_AthaGene_InOrder");
//		List<List<String>> map = cmmFun.readMapFile("/home/setas/git/IndependentParsers/Independent Parsers/src/files/GF_GO_Desc_allHierarchy_28Eudicot");
//		
//		
//		for(String gfID : EudGFsInOrder){			
//			
//			for(List<String> ls : map){
//				
//				if(ls.get(0).equals(gfID)){
//					
//					System.out.println(ls.get(0)+"\t"+ls.get(1)+"\t"+ls.get(2));
//				}
//				
//			}
//			
//		}
		
		ArrayList<String> GFsInOrder_37spe = cmmFun.read1ColFile_String("/home/setas/git/IndependentParsers/Independent Parsers/src/files/9178coreGF_inOrderLam");
		List<List<String>> map = cmmFun.readMapFile("/home/setas/git/IndependentParsers/Independent Parsers/src/files/genefamily.angiosperm.dup.rootsupport70.paml.txt");
		
		for(String gfID : GFsInOrder_37spe){			
			
		for(List<String> ls : map){
			
//			System.out.println("gfID: "+gfID +"  ls(0): "+ ls.get(0));
			
			if(ls.get(0).equalsIgnoreCase(gfID)){
				
//				System.out.println("gfID: "+gfID +"  ls(0): "+ ls.get(0));
				
				for(int i=0; i< ls.size(); i++){
				System.out.print(ls.get(i)+"\t");}
				
				System.out.println();
			}
			
		}
		
	}
		
	}
}
