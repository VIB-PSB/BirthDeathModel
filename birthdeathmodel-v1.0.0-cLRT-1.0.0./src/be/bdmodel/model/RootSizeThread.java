package be.ugent.psb.setas.bdmodel.model;

import java.io.IOException;
import java.util.Properties;


public class RootSizeThread extends Thread {
		
	private CuttingPlaneMethod cpm;
	

	public RootSizeThread(CuttingPlaneMethod cpm) {
		this.cpm = cpm;
	}

	@Override
	public synchronized void run() {
		
		cpm.findOptimalLambda();
		
//		return optimalLambda; should be able to return [optimal lambda, optimal lk] such that later we can loop over this output and pick the maximum lk 
		// and corresponding lambda and root size for one family

		System.out.print(cpm.getRootSize()+"\t"+cpm.getOptimalLambda()+"\t"+cpm.getFoptimalLambda());
		System.out.println();

	}
	
	
}
