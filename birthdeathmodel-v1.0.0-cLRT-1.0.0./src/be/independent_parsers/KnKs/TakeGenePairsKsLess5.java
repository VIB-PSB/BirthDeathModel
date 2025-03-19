package be.ugent.psb.setas.independent_parsers.KnKs;

import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class TakeGenePairsKsLess5 {


	public static void main(String[] args) {

		CommonFunctions cmmFunct = new CommonFunctions();

		List<List<String>> submap = cmmFunct
				.subMap("/home/setas/git/IndependentParsers/Independent Parsers/src/files/genefamily.angiosperm.dup.rootsupport70.paml.Atha.txt",
						9, 5.0);
		
		for(List<String> ls: submap){
			
			for(int i=0; i<ls.size(); i++){
				
			
			if(i!= 2 && i != 3){	
			System.out.print(ls.get(i)+"\t");
				}
			
			}
			System.out.println();
		}

	}

}
