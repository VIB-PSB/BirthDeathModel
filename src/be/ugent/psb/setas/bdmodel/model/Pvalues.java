
package be.ugent.psb.setas.bdmodel.model;

/*
 * #%L
 * BirthDeathModel
 * %%
 * Copyright (C) 2017 VIB/PSB/UGent - Setareh Tasdighian
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

import java.util.ArrayList;

/**
 * p-values are the number of generated Observations that have less
 * likelihood than the reference likelihood + half of the ones with equal
 * likelihood divided by the whole number of generated observations
 */
public class Pvalues {


	private Node root;
	private ArrayList<Node> arraylistOfNodes;
	private int numberOfObservations;
	private GenerateRandomGeneCountProfiles generateObservations;
	private CalculateLikeLihoods likelihood;
	private double lambda;

	static final double EPSILON = 1e-60;

	public Pvalues(Node root, ArrayList<Node> arrayylistOfNodes , double lambda, int numberOfRandomGeneratedGeneFamilyProfiles){
		
		System.out.println("Warning: Class Pvalues is not using any Cache");
		this.root = root;
		this.arraylistOfNodes = arrayylistOfNodes;
		this.numberOfObservations = numberOfRandomGeneratedGeneFamilyProfiles;
		this.generateObservations = new GenerateRandomGeneCountProfiles(0,this.root.getMaxGeneCountAtNode(),false);
		this.likelihood = new CalculateLikeLihoods(lambda,arrayylistOfNodes.get(0).getMaxGeneCountAtNode()+1);	
	}
	
	public Pvalues(Node root, ArrayList<Node> arrayyListOfNodes , double lambda, int numberOfObservations, TransitionProbabilityCalculator probabilityCalculator, int lengthOfMCMC ){
		this.root = root;
		this.arraylistOfNodes = arrayyListOfNodes;
		this.lambda= lambda;
		this.numberOfObservations = numberOfObservations;
		this.generateObservations = new GenerateRandomGeneCountProfiles(5,root.getMaxGeneCountAtNode(),false,probabilityCalculator,lengthOfMCMC);
		this.likelihood = new CalculateLikeLihoods(lambda,arrayyListOfNodes.get(0).getMaxGeneCountAtNode()+1, probabilityCalculator);	
	}
	

	public static boolean isEqual(double double1, double double2, double precision){
		return (Math.abs(double1-double2) < precision);
	}
	
	
	public double calculateConditionalPvalues(int geneCountAtNode, int rootSize, double referenceLogLikelihood) {
		
		    double pValue =0;
		
			for (int randomObservation = 1; randomObservation <= numberOfObservations; randomObservation++) {
				
				int[] generatedObservationArray = generateObservations.generateGeneCountProfile(root, rootSize, this.lambda);
				
				root.setLeafValues(generatedObservationArray);

				double[] likelihoodOfThisObservation = likelihood.calculateInternalLikelihoods(arraylistOfNodes);	
				
				double[] logLikelihoodOfThisObservation = MathematicalOperations.giveLogarithm10Array(likelihoodOfThisObservation);
				
				if(isEqual(logLikelihoodOfThisObservation[geneCountAtNode],referenceLogLikelihood,EPSILON)){
					
					pValue += 0.5;			
				}
						
				else if (logLikelihoodOfThisObservation[geneCountAtNode] < referenceLogLikelihood) {
					pValue += 1;		
				}
				
			}
			pValue = (pValue*1.0)/(numberOfObservations*1.0);
		
		return pValue;
	}

}
