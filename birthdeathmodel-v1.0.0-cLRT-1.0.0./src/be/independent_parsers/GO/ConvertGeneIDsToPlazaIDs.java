package be.ugent.psb.setas.independent_parsers.GO;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class ConvertGeneIDsToPlazaIDs {
	
//	public List<String> readSpeciesTagsToConvertToPlazaIDs(String filename){
//		
//		FileReader fr = null;
//		try {
//			fr = new FileReader(filename);
//		} catch (FileNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//
//		Scanner sc = new Scanner(fr);
//		List<String> table = new ArrayList<String>();
//
//		while (sc.hasNextLine()) {
//			
//			String line = sc.nextLine().trim();
//			table.add(line);
//		}
//		sc.close();
//		
//		return table;
//	}
	
	public static void main(String [] args){
		
		
		CommonFunctions cmmf = new CommonFunctions();
//		ConvertGeneIDsToPlazaIDs crvtGIDs = new ConvertGeneIDsToPlazaIDs();
		
//		String file1 = "/home/setas/Desktop/AllGenes_28spePlazaID";
//		String file1 = "/home/setas/Desktop/AllGenes_AllSpecies_AllGFsNEW";
		String file1 = "/home/setas/Desktop/AllGenes_AllGFs_36spePlazaIDs_Tcac";
		
		ArrayList<String> myGFID = cmmf.readColX_String(file1,0);
		ArrayList<String> myGeneID = cmmf.readColX_String(file1,1);
		
		String file2="/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/TBpercentages/conversion_Plaza/gra.filtered.tab"; //Conversion File From Plaza		
//		ArrayList<String> PlazaGeneID = cmmf.readColX_String_SemiColon(file2,0);
//		ArrayList<String> nonPlazaID = cmmf.readColX_String_SemiColon(file2,1);
//		ArrayList<String> nonPlazaID = cmmf.readColX_String_SemiColon(file2,2); // only in case of Sotub and evm Hvul
//		ArrayList<String> nonPlazaID = cmmf.readColX_String_SemiColon(file2,3); //only for macu
		String speciesTag = "Cotton"; 
		ArrayList<String> PlazaGeneID = cmmf.readColX_String(file2,1); //only for tcac and grai
		ArrayList<String> nonPlazaID = cmmf.readColX_String(file2,0);
		
//		List<String> myTags = crvtGIDs.readSpeciesTagsToConvertToPlazaIDs("");
						
		
		for(int i=0; i<myGeneID.size(); i++){
			
			
			if(myGeneID.get(i).split(speciesTag).length > 1) {
				
				/** only for csat **/
//				String [] parts = myGeneID.get(i).split("M");
//				String newGeneID =(parts[0]+"G"+parts[1]);				
//				System.out.println(myGFID.get(i)+"\t"+ newGeneID);
				
//				String secondPartOfID = myGeneID.get(i).split(speciesTag)[1]; //only for Cmel		
//				int index = cmmf.searchListString_index("MELO"+secondPartOfID, nonPlazaID);				
//				int index = cmmf.searchListString_index(myGeneID.get(i)+".g", nonPlazaID);
				
				int index = cmmf.searchListString_index(myGeneID.get(i), nonPlazaID);
					
					if(index >= 0){
					String plazaID = PlazaGeneID.get(index);		
					System.out.println(myGFID.get(i)+"\t"+plazaID);}
					
					else{
						System.out.println(myGFID.get(i)+"\t"+"No Plaza ID Found: "+myGeneID.get(i));
					}
			}
			
			else{
				System.out.println(myGFID.get(i)+"\t"+myGeneID.get(i));
			}
		}
		
	}

}
