package be.ugent.psb.setas.independent_parsers.KnKs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

import be.ugent.psb.setas.bdmodel.model.MathOperations;
import be.ugent.psb.setas.independent_parsers.CommonFunctions;

public class CreatKsBinsForKn {

	public static void main(String[] args) {

		CommonFunctions cmm = new CommonFunctions();

		String path = 
        "/mnt/shares/biocomp/groups/group_esb/setas/Project1/Results/9178coreGF-CPMpv-tier1-37speMGCF5/KnKs/GFsOrder_genDup/Ks5_minKs/TopBottom/GenPairs_PAML/Zmay_Ks5_Top709.txt";
		
		ArrayList<List<Double>> allBins = new ArrayList<List<Double>>();
		double ksMax =3.00;
		double binSize =0.1;
		
		int numOfBins = (int)(Math.floor(ksMax/binSize)) ;
		
		for(int i=0; i<numOfBins; i++){
			
			List<Double> lsd = new ArrayList<Double>();
			allBins.add(lsd);
		}

		ArrayList<String> gfId = cmm.readColX_String(path, 0);
		ArrayList<String> gene1 = cmm.readColX_String(path, 1);
		ArrayList<String> gene2 = cmm.readColX_String(path, 2);
		ArrayList<Double> t = cmm.readColX_double(path, 3);
		ArrayList<Double> s = cmm.readColX_double(path, 4);
		ArrayList<Double> n = cmm.readColX_double(path, 5);
		ArrayList<Double> omega = cmm.readColX_double(path, 6);
		ArrayList<Double> Kn = cmm.readColX_double(path, 7);
		ArrayList<Double> Ks = cmm.readColX_double(path, 8);

		for (int i = 0; i < Ks.size(); i++) {

			if(Ks.get(i)<= ksMax){
			int bin = (int) (Math.floor(Ks.get(i) / binSize));
			
			allBins.get(bin).add(Kn.get(i));
			}

//			System.out.println(gfId.get(i)+"\t"+gene1.get(i)+"\t"+gene2.get(i)+"\t"+t.get(i)+"\t"+s.get(i)+"\t"+n.get(i)+"\t"+omega.get(i)+"\t"+Kn.get(i)+"\t"+Ks.get(i) + "\t" + bin);
		}
		
		double [] means =new double[allBins.size()];
		double [] stds = new double [allBins.size()];
		
		for(int i=0; i<allBins.size(); i++){
			
			List<Double> oneBin = allBins.get(i);
			
			double mean = MathOperations.calcAverage(oneBin);
			double std = MathOperations.calcStd(oneBin);
			
			means[i]= mean;
			stds[i] = std;
			
		}

		for(int k=0; k<means.length;k++){
			System.out.println(means[k]+"\t"+stds[k]);
		}
		
	}

}
