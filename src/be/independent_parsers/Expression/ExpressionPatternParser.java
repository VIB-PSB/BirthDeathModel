package be.ugent.psb.setas.independent_parsers.Expression;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

public class ExpressionPatternParser {
	
	private ArrayList<ArrayList<String>> gfEntries = new ArrayList<ArrayList<String>>();
	private Map<Couple, ArrayList<String>> coupleMap = new HashMap<Couple, ArrayList<String>>();
	
	private class Couple{		
		private HashSet<String> probeSet = new HashSet<String>();
		
		public Couple(String p1, String p2){
			probeSet.add(p1);
			probeSet.add(p2);
		}

		@Override
		public int hashCode() {
			return probeSet.hashCode();
		}

		@Override
		public boolean equals(Object obj) {
			if (!(obj instanceof Couple)) return false;
			Couple couple = (Couple)obj;
			
			return this.probeSet.equals(couple.probeSet);
		}
				
	}
	
	private void parseProbeList(File geneFamilyFile) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(geneFamilyFile));
		String line = in.readLine();
		while (line != null){
			String[] cols = line.split("\t");
			String gf = cols[0];
			String probe1 = cols[2].split("\\|")[1];
			String probe2 = cols[3].split("\\|")[1];
			String ks = cols[9];
			
			ArrayList<String> entry = new ArrayList<String>();
			entry.add(gf);
			entry.add(probe1);
			entry.add(probe2);
			entry.add(ks);
			
			gfEntries.add(entry);
			coupleMap.put(new Couple(probe1, probe2), entry);
			
			line = in.readLine();
		}
		in.close();
	}
	
	private void parseCornetFile(File cornetFile) throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(cornetFile));
		String line = in.readLine();
		while (line != null){
			String[] cols = line.split("\t");
			String probe1 = cols[0];
			String probe2 = cols[1];
			String exprDiv = cols[2];
			
			Couple key = new Couple(probe1, probe2);
			
			ArrayList<String> entry = coupleMap.get(key);
			if (entry != null){
				entry.add(exprDiv);
			}
			
			line = in.readLine();
		}
		in.close();
	}
	
	private void writeOutput(File outputFile) throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		for (ArrayList<String> entry : gfEntries){
			if (entry.size() == 4){
				entry.add(new String("N/A"));
			}
			for (int i = 0; i < entry.size(); i++){
				String col = entry.get(i);
				out.write(col);
				if (i<entry.size()-1){
					out.write("\t");					
				} else {
					out.newLine();
				}
			}
		}
		out.close();
	}
	
	public static void main(String[] args) {
		File gfFile = new File("/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/newfromLor/genefamily.angiosperm.core.txt.codeml");
		File cornetFile = new File("/home/setas/Desktop/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/ExpPatterns/newfromLor/CornetBYadj.txt");
		
		File outputFile = new File("/home/setas/Desktop/cornetscan.txt");
		
		ExpressionPatternParser parser = new ExpressionPatternParser();
		
		try {
			parser.parseProbeList(gfFile);
			parser.parseCornetFile(cornetFile);
			parser.writeOutput(outputFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

}
