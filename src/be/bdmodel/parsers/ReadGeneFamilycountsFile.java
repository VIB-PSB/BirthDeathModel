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
 * Reads the gene counts file which is a tab delimited file
 * with headers as gene family IDs, number of genes, number of species
 * and then all the gene counts of all species in the phylogeny
 * 
 * CAUTION: The columns of this file must be in the same order as the species appear in the newick format tree file.
 * CAUTION: The columns of the tab separated gene count file are as follows:
 * Column 0: GeneFamilyID , Column 1: number of genes, Column 2: number of Species, Coulmns :3,4,...,40 : gene counts of all species
 */
public class ReadGeneFamilycountsFile {

	private ArrayList<String> geneFamilyIDs;

	public ArrayList<String> getGfIDs() {
		return geneFamilyIDs;
	}
	
	public List<List<Integer>> readGeneCountsFile(String geneCountsFile) {

		FileReader fileReader = null;
		try {
			fileReader = new FileReader(geneCountsFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		Scanner scanner = new Scanner(fileReader);
		scanner.nextLine();

		List<List<Integer>> geneFamilyCounts = new ArrayList<List<Integer>>();

		geneFamilyIDs = new ArrayList<String>();

		while (scanner.hasNextLine()) {

			String line = scanner.nextLine();
			String[] parts = line.split("\t");

			List<Integer> geneCounts = new ArrayList<Integer>();
			for (int i = 3; i < parts.length; i++) {
				int geneCountsValue = Integer.parseInt(parts[i]);
				geneCounts.add(geneCountsValue);
			}
			geneFamilyCounts.add(geneCounts);
			geneFamilyIDs.add(parts[0]);
		}
		scanner.close();
		return geneFamilyCounts;
	}

}
