package be.ugent.psb.setas.bdmodel.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Generates subsets of given size from an initial starting set
 * @param <T> Type of List items
 */

public class ChooseKoutOfN<T> {
	private List<T> inputSet;
	List<List<T>> outputList = new ArrayList<List<T>>();

	/**
	 * Construct subset generator with initial set
	 * 
	 * @param inputList
	 */
	public ChooseKoutOfN(List<T> inputList) {
		this.inputSet = new ArrayList<T>(inputList);
	}

	/**
	 * Generate all subsets of inputSet with size k
	 *
	 * @param k
	 *            size of the subsets
	 * @return list of all possible subsets
	 */
	public List<List<T>> generateSubSets(int k) {
		subset(k, 0, new ArrayList<T>());
		return outputList;
	}

	private void subset(int left, int index, ArrayList<T> l) {
		if (left == 0) {
			this.outputList.add(new ArrayList<T>(l));
			return;
		}
		for (int i = index; i < this.inputSet.size(); i++) {
			l.add(inputSet.get(i));
			subset(left - 1, i + 1, l);
			l.remove(l.size() - 1);
		}
	}

	/**
	 * Write results to standard out
	 */
	public void writeOutput() {
		for (List<T> mySet : this.outputList) {
			for (T item : mySet) {
				System.out.print(item + "  ");
			}
			System.out.println();
		}

	}

}
