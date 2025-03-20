package be.ugent.psb.setas.bdmodel.parsers;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ugent.psb.setas.bdmodel.model.Node;

public class BuildTree {

	/**
	 * NOTE: for type I input, be.ugent.psb.setas.bdmodel.test data, e.g: R,A,5 The Newick format trees
	 * are built through Newick parser read the list of nodes and branches
	 * between them from input file and builds the tree
	 * 
	 * @param filename
	 * @return
	 */

	public List<Node> readInput(String filename) {

		FileReader fr = null;
		try {
			fr = new FileReader(filename);
		} catch (FileNotFoundException e) {e.printStackTrace();}

		Scanner sc = new Scanner(fr);

		/** table will contain all built nodes from input lines **/
		List<Node> table = new ArrayList<Node>();
		Node node = new Node();
		Node root = new Node();

		while (sc.hasNextLine()) {
			String line = sc.nextLine();
			String[] names = line.split(",");

			/** make a node and its children from each line of input: **/
			root = node.AddChildren(names[0], names[1], names[2],
					Double.parseDouble(names[3]), Double.parseDouble(names[4]));

			table.add(root);
		}
		sc.close();
		return table;
	}

	public Node searchList(List<Node> l, String name) {
		int size = l.size();
		for (int i = 0; i < size; i++) {
			if (l.get(i).getName().equals(name)) {return l.get(i);}
		}
		return null;
	}

	/** the table which is built in readInput should be passed to ln: */
	public Node buildTree(List<Node> ln, Node currentRoot) {
		
		Node tree = currentRoot;

		/** to see if the new input line- new Node- does not already exist, and add its
		 children instead of creting a separate node **/
		Node left = searchList(ln, tree.getLeftChild().getName());
		Node right = searchList(ln, tree.getRightChild().getName());

		// since leaves will produce exception:
		if (left != null) {
			double leftblen = tree.getLeftChild().getbLen();
			tree.setLeftChild(buildTree(ln, left));
			tree.getLeftChild().setbLen(leftblen);
			
		} else {
			tree.getLeftChild().isLeaf = true;}

		if (right != null) {
			
			double rightblen = tree.getRightChild().getbLen();
			tree.setRightChild(buildTree(ln, right));
			tree.getRightChild().setbLen(rightblen);
		}

		else {tree.getRightChild().isLeaf = true;}

		return tree;
	}

//	public static void main(String args[]) {
//		BuildTree bt = new BuildTree();
//		List<Node> list = bt
//				.readInput("/home/setas/workspace/BirthDeathModel/src/Files/plaza.txt");
//		// System.out.println(n.size());
//		Node tree = bt.buildTree(list, list.get(0));
//
//		System.out.println(tree.getName()); // root
//		System.out.println(tree.getLeftChild()); // A
//		System.out.println(tree.getLeftChild().getbLen()); // 10
//
//		System.out.println(tree.getLeftChild().getLeftChild().getbLen()); // 8
//		System.out.println(tree.getLeftChild().getLeftChild().isLeaf); // true
//		System.out.println(tree.getLeftChild().isLeaf); // false
//	}

}
