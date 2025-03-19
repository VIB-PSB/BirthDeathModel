package be.ugent.psb.setas.independent_parsers.RVSEng;

import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class CalcStatRvsEng_ParsePythonNiFiles {
	
	
	public static boolean isParsable(String input){
	    boolean parsable = true;
	    try{
	        Integer.parseInt(input);
	    }catch(NumberFormatException e){
	        parsable = false;
	    }
	    return parsable;
	}
	
	public static boolean isParsable2(String input){
	    boolean parsable = true;
	    try{
	        Double.parseDouble(input);
	    }catch(NumberFormatException e){
	        parsable = false;
	    }
	    return parsable;
	}
	
	public static void main(String[] args) {

		CommonFunctions cmf = new CommonFunctions();
		
		String path_python_ni_file = "/home/setas/Desktop/setas/ReverseEngineering/Musa3TauMon2/H0isRmWGD/pythonNi_Conventional_H0RMWGD";
		
		List<List<String>> map_GFid_Ni = cmf.readMapFile(path_python_ni_file);
		
		for(List<String> ls :map_GFid_Ni){
			
			System.out.print(ls.get(0)+"\t");			
			
			for(int i=2; i<ls.size();i++){ //i=1 = rootSize
				
				String ni= ls.get(i);

				if(isParsable2(ni)){
				
				double ni_double = Double.parseDouble(ni);
				double ni_LessThanQstar = ni_double;
				
				if(ni_LessThanQstar > 900){
			
					int numOfWGD = i-2;
					System.out.print(numOfWGD+ "\t");
				}
				
				}
			}
			
			System.out.print("\n");
		}
		
		
	}
	
	

}
