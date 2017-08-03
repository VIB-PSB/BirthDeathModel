package be.ugent.psb.setas.bdmodel.parsers;

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

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

/**
 * reads a file containing all WGDs on the tree
 * 
 * @param filename
 * @return a list of WGDs to be built via addWGD and InsertWGD in class Node.
 * 
 *         Warning: If there are multiple WGMs on a branch, the order in which
 *         they appear in the text file is important: first the older events and
 *         then the younger ones
 */
public class WGMparser {

	/**
	 * Read the input file for WGMs:
	 * every line should look like: WGD/T,NameOfParentNode, NameOfChildNode, branchLength
	 * @param wgmFile
	 * @return
	 */
	public List<List<String>> readWGMfile(String wgmFile) {

		FileReader fileReader = null;
		try {
			fileReader = new FileReader(wgmFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner scanner = new Scanner(fileReader);

		List<List<String>> listOfAllWGMs = new ArrayList<List<String>>();

		while (scanner.hasNextLine()) {

			String line = scanner.nextLine();
			String[] parts = line.split(",");

			List<String> oneWGM_parentName_childName_BranchLength = new ArrayList<String>();

			for (int i = 0; i < parts.length; i++) {
				oneWGM_parentName_childName_BranchLength.add(parts[i]);
			}
			listOfAllWGMs.add(oneWGM_parentName_childName_BranchLength);

		}
		scanner.close();
		return listOfAllWGMs;
	}

}
