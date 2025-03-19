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

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class CommonFunctions {

	/**
	 * Reads a tab delimited file and returns the file as a list of list of integers, each list of integer being one line.
	 * @param tabDelimitedFile
	 * @return
	 */
	public List<List<String>> readTabDelimitedFile(String tabDelimitedFile) {

		FileReader fileReader = null;
		try {
			fileReader = new FileReader(tabDelimitedFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner scanner = new Scanner(fileReader);
		List<List<String>> allLines = new ArrayList<List<String>>();

		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();

			if (!line.isEmpty()) {

				String[] parts = line.split("\t");
				List<String> listOfStringsInLine = new ArrayList<String>();

				for (int i = 0; i < parts.length; i++) {
					listOfStringsInLine.add(parts[i]);
				}

				allLines.add(listOfStringsInLine);
			}
		}
		scanner.close();
		return allLines;
	}
	

}
