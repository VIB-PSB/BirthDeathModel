package be.ugent.psb.setas.bdmodel.model.singlecopyproject;

import be.ugent.psb.setas.bdmodel.model.Node;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;

public class GenerateNullDistributionLambda {

	public void generateRandomGeneCountFile(String treeFile, String WGDFile, double partitionSize,
			int defaultmaxNodeSize, int numOfObservations) {

		Node root = SpeciesTreeParser.buildInsertWGDsandPartitionTree(treeFile, WGDFile, partitionSize,
				defaultmaxNodeSize);
	}
}
