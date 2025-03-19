package be.ugent.psb.setas.bdmodel.model;

import java.util.List;

import be.ugent.psb.setas.bdmodel.model.RVS_Engineering.MGFsimulationsRemovingEachWGD;

public class ClusterGFrecords {
	
	private double partitionSizeLambda;
	private double maxLambda;
	
	public ClusterGFrecords(double partitionSizelambda, double maxLam){
		
		this.partitionSizeLambda = partitionSizeLambda;
		this.maxLambda = maxLam;
	}

	public double [] calcIntervals(int index){
		
		double [] interval = new double[2];
		
		interval[0] = (index-1)*0.01;
		interval[1] = index *0.01;
		
		return interval;
	}
	
	public static boolean isEqual(double a, double b, double PRECISION){
		return (Math.abs(a-b) < PRECISION);
	}
	
	public boolean ifLambdaInInterval(double [] interval, double lam){
		
		boolean b= false;
		
		if ((lam > interval[0] || isEqual(lam,interval[0], 1e-6) )&& lam < interval[1]){
			
			b= true;
		}
		
	
		return b;
	}
	
	public static void main(String[] args) {

		ProbCalculator probCalc = new ProbCalculator();
		MGFsimulationsRemovingEachWGD mgfsim = new MGFsimulationsRemovingEachWGD(
				probCalc);
//		
//		double partitionSizeLambda = 0.01;
//		double maxLam = 3;
//		
//		int numberOfIntervals = (int) Math.floor(maxLam /partitionSizeLambda);	
//		ClusterGFrecords clsGF = new ClusterGFrecords(partitionSizeLambda,maxLam);
//			
//		for(int i=1; i<numberOfIntervals;i++){
//			
//			double [] interval = clsGF.calcIntervals(i);		
//		}
//
		String pathLRTfile = "/home/setas/Desktop/sorted_TauMon2Musa3_Rm/sorted_TauMon2Musa3_RmWGDs.txt";
		List<List<Double>> LRTtable = mgfsim.readMGFfile(pathLRTfile); 
//
		for (int mgf = 0; mgf < 9178; mgf++) {

			String gfID = mgfsim.getGfIDs().get(mgf);

			List<Double> mgfRecord = LRTtable.get(mgf);
			
			/** For lambda **/
//			double lambdaOriginal = mgfRecord.get(1);
//			int numberOfInterval = (int)(lambdaOriginal/partitionSizeLambda);
//			System.out.println(numberOfInterval);

			
			
		}
		
		
		
	}
}
