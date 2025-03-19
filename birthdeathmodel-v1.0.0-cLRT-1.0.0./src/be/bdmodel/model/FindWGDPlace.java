package be.ugent.psb.setas.bdmodel.model;

import java.util.ArrayList;

import be.ugent.psb.setas.bdmodel.parsers.NewickParser;
import be.ugent.psb.setas.bdmodel.parsers.SpeciesTreeParser;

public class FindWGDPlace {

	public Node root;
	
	public Node findNodeWithName(String name){
		
		ArrayList<Node> aryln = root.postOrder(root);
		Node found = new Node();
		
		for (int i=0; i<aryln.size();i++){
			
			if(aryln.get(i).getName().equals(name)){
			
				found = aryln.get(i);
			}
			
		}
		
		if(found.getName() == null){
			
			System.out.println("A leaf with name "+name+ " does not exist in this tree. ");
		}
		
		return found;
	}
	
	public void placeWGMonBranch(int multFactor, double distanceWGDToLeaf, Node leaf){

		double distanceLeafToRoot= leaf.calcDistToRoot();
		
		if(distanceWGDToLeaf > distanceLeafToRoot){
			
			System.out.println("This distance of WGD to leaf node is more than the distance to root!");
		}

		double sumBlens = 0;
		Node Child = leaf;
		Node Parent = Child.getParent();
		
		sumBlens += Child.getbLen();
		
//		System.out.println("sumBlen "+sumBlens);
		
		if(sumBlens < distanceWGDToLeaf){
						
			System.out.println("distance WGD to leaf: "+ distanceWGDToLeaf + "  sumBlens: "+ sumBlens);
			placeWGMonBranch(multFactor, (distanceWGDToLeaf-sumBlens),Parent);
			
			
		}
		
		else{
			
			double WGMBlen = sumBlens - distanceWGDToLeaf;
			
//			System.out.println("WGM Blen: "+WGMBlen);
			
			if(multFactor==2){
			String S = "WGD,"+ Parent.getName()+","+Child.getName()+","+WGMBlen;
			System.out.println(S);
			
			double S2 = WGMBlen/ Child.getbLen() ;
			System.out.println("index: "+S2);
			System.out.println();
			}
			
			if(multFactor==3){
				String S = "WGT,"+ Parent.getName()+","+Child.getName()+","+WGMBlen;
				System.out.println(S);
			}		
		}	
	}
	
	public static void main(String[] args){
		
//		NewickParser np = new NewickParser();
//		Node root = np.buildAndPartitionTree("/home/setas/workspace/BirthDeathModel/src/Files/Trees/37spe-MGCF5.txt",100);
		
		double partitionSize = 0;
		int defaultmaxNodeSize =100;
		Node root = SpeciesTreeParser.buildAndPartitionTree_NoWGD(args[0], partitionSize, defaultmaxNodeSize);
		
		FindWGDPlace fwp = new FindWGDPlace();
		fwp.root = root;
		
		Node leaf = fwp.findNodeWithName("Ljap");
		fwp.placeWGMonBranch(2,0.7566,leaf);
		
	}
}
	