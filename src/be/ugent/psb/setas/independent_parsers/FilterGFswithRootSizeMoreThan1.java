package be.ugent.psb.setas.independent_parsers;

import java.util.List;

public class FilterGFswithRootSizeMoreThan1 {
	
	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();

		String path ="/home/setas/Desktop/MGFfileForGFsMore10wgdRvs";
		
		List<List<String>> map= cmf.readMapFile(path);
		
        for (int i=0; i< map.size(); i++){
        
        	List<String> currentRecord = map.get(i);
        	
        	if((int)(Double.parseDouble(currentRecord.get(1))) == 1){
        		
        		for(int j=0; j<currentRecord.size(); j++){
        		System.out.print(currentRecord.get(j)+"\t");
        		}
        		System.out.println();
        	}
        	
        }
        }
		
}
