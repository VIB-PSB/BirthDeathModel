package be.ugent.psb.setas.bdmodel.parsers;

import java.io.FileReader;
import java.util.ArrayList;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Scanner;
import java.util.Stack;

import be.ugent.psb.setas.bdmodel.model.Node;

public class NewickParser {

	public String[] convertQofStringToArray(Queue<String> q) {
		int size = q.size();
		String[] b = new String[size];

		for (int k = 0; k < size; k++) {
			b[k] = q.remove();
		}
		return b;
	}

	// without paranthesis: (also without , and : and ; )
	public String[] readInput_removeChar_InputString(String filename) {

		FileReader fr = null;
		try {
			fr = new FileReader(filename);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner sc = new Scanner(fr);
		String line = new String();

		Queue<String> q = new LinkedList<String>();

		while (sc.hasNextLine()) {
			line = sc.nextLine();
			String[] s = line.split("(\\,)|(:)|(\\()|(\\))|(;)");

			for (int i = 0; i < s.length; i++) {
				q.add(s[i]);
			}
		}
		sc.close();

		String[] nodeInfo = convertQofStringToArray(q);

		return nodeInfo;
	}

	// feed output of readInput
	public ArrayList<Node> makeNodes_setName_blen(String[] nInfo) {

		int len = nInfo.length;

		Queue<String> qs = new LinkedList<String>();
		ArrayList<Node> arn = new ArrayList<Node>();

		for (int k = 0; k < len; k++) {
			if (!nInfo[k].isEmpty()) {
				qs.add(nInfo[k]);
			}
		}

		String[] nodeInfoEdit = convertQofStringToArray(qs);

		int lent = nodeInfoEdit.length;
		for (int i = 0; i < lent - 1; i = i + 2) {
			Node n = new Node();
			n.setName(nodeInfoEdit[i]);
			n.setbLen(Double.parseDouble(nodeInfoEdit[i + 1]));
			arn.add(n);
					}
		Node r = new Node();
		r.setName(nodeInfoEdit[lent - 1]);
		arn.add(r);
		return arn;
	}

	public static boolean isNumeric(String str) {
		return str.matches("-?\\d+(\\.\\d+)?"); // match a number with optional
												// '-' and decimal.
	}

	// WITH paranthesis, without numbers: (without , and : and ; )
	public String[] readInput_KeepParanthesis_rmRest(String filename) {

		FileReader f = null;
		try {
			f = new FileReader(filename);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Scanner scan = new Scanner(f);
		String line = new String();

		Queue<String> q = new LinkedList<String>();

		while (scan.hasNextLine()) {
			line = scan.nextLine();
			String[] s = line.split("(\\,)|(:)|(;)");

			for (int i = 0; i < s.length; i++) {
				if (!isNumeric(s[i])) {
					q.add(s[i]);
				}
			}
		}
		scan.close();

		String[] nodeInfoTwo = convertQofStringToArray(q);

		return nodeInfoTwo;
	}

	// feed output of readInputTwo to this
	public int[] calculateDepth(String[] treeString2) {

		int lenx = treeString2.length;
		int[] b = new int[lenx];

		int d = 0;

		for (int i = 0; i < lenx; i++) {
			String s = treeString2[i];
			
			int lenOfSlot = s.length();
			char[] ch = new char[lenOfSlot];

			for (int j = 0; j < lenOfSlot; j++) {

				ch[j] = s.charAt(j);
				// System.out.println("chars: "+ i+""+ j+""+ ch[j]);
				if (ch[j] == '(') {
					d += 1;
				}
				if (ch[j] == ')') {
					d -= 1;
				}
			}
			b[i] = d;
		}
		 for(int i=0;i<b.length;i++){
		 }
		return b;
	}

	public void setDepth(ArrayList<Node> arrayListNodes, String[] treeString2) {

		int[] depths = calculateDepth(treeString2);
		int size = arrayListNodes.size();
		
		for (int i = 0; i < size; i++) {
			
			arrayListNodes.get(i).depth = depths[i];			
		}
	}

	// w.r.t the structure of newick format: first node: left most node and so
	// on
	public Node setParentChildRelations(ArrayList<Node> arrayListNodes) {
		Node root = new Node();
		Stack<Node> st = new Stack<Node>();

		st.push(arrayListNodes.get(0));

		for (int i = 1; i < arrayListNodes.size(); i++) {

			Node current = arrayListNodes.get(i);

			if (current.depth >= st.peek().depth) {
				st.push(current);
			}

			else {
				st.peek().setParent(current);
				current.setRightChild(st.pop());
				
				st.peek().setParent(current);
				current.setLeftChild(st.pop());
				st.push(current);
			}

		}

		root = st.pop();
		root.isRoot =true;
		return root;
	}
	
	public void setDistanceToRoot(ArrayList<Node> arrayListNodes) {

		int size = arrayListNodes.size();
		
		for (int i = 0; i < size; i++) {

			if(arrayListNodes.get(i).isRoot == false){
			
			double distanceToRoot= arrayListNodes.get(i).calcDistToRoot();
			
			arrayListNodes.get(i).setDistToRoot(distanceToRoot);
			}
			
			else{
				arrayListNodes.get(i).setDistToRoot(0);
			}
		}
	}

	public void setIsleafBoolean(Node myNode) {

		if (myNode != null) {
			Node left = myNode.getLeftChild();
			Node right = myNode.getRightChild();

			if (right == null && left == null) {
				myNode.isLeaf = true;
			}

			setIsleafBoolean(left);
			setIsleafBoolean(right);
		}
	}
	
	public void setDefaultMaxNodeSize(ArrayList<Node> arrayListNodes, int defaultMaxNodeSize) {
		for (int i = 0; i < arrayListNodes.size(); i++) {
			arrayListNodes.get(i).setmaxNodeSize(defaultMaxNodeSize);
		}
	}
	
	public Node buildTree(String newickTreefilename , int defaultMaxNodeSize) {

		String[] treeString1 = readInput_removeChar_InputString(newickTreefilename);
		String[] treeString2 = readInput_KeepParanthesis_rmRest(newickTreefilename);

		ArrayList<Node> arrayListNodes = makeNodes_setName_blen(treeString1);
		setDepth(arrayListNodes, treeString2);
	
		Node r = setParentChildRelations(arrayListNodes);
		setDistanceToRoot(arrayListNodes);
			
		setDefaultMaxNodeSize(arrayListNodes,defaultMaxNodeSize);
		
		setIsleafBoolean(r);
		
//		ArrayList<Node> aryl = r.getLeaves();
//		System.out.println(aryl.size());
//		for(Node l : aryl){System.out.println(l.getName());}
		
//		ArrayList<Node> allNodes = root.postOrder(root);
//		for(Node n :allNodes){System.out.println(n.getName());}
		
		return r;
	}

//	 public static void main(String[] args) throws IOException {
//	
//	 NewickParser np = new NewickParser();
//	 Node root = np
//	 .buildTree("/home/setas/workspace/BirthDeathModel/src/Files/toyTree.txt");
//	
//	 WGM wgd = new WGM();
//	 List<List<String>> wgdList =
//	 wgd.readInputFile("/home/setas/workspace/BirthDeathModel/src/Files/wgd.txt");
//	 root.insertWGM(root, wgdList);
//	
//	 ArrayList<Node> al = root.getLeaves();
//	
//	 root.setNumberOfLeaves(al.size());
//	
//	 ArrayList<Node> a= root.postOrder(root);
//	 
//	 for(int i=0; i<a.size();i++){
//		 
//		 System.out.println(a.get(i).getName()+"  "+ a.get(i).getbLen());
//		 
//	 }
//	
//	 }

}
