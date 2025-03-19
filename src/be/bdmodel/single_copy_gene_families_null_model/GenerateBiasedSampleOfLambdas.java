package be.ugent.psb.setas.bdmodel.model.singlecopyproject;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;
import java.util.Scanner;

public class GenerateBiasedSampleOfLambdas {

	int sampleSize;

	public int[] readFrequencies(String fileName) {

		FileReader fin = null;
		try {
			fin = new FileReader(fileName);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		/* The first line is a header: */
		Scanner sc = new Scanner(fin);
		List<Integer> table = new ArrayList<Integer>();

		while (sc.hasNextLine()) {
			String line = sc.nextLine();

			if (!line.equalsIgnoreCase("")) {

				int parsed = Integer.parseInt(line);

				table.add(parsed);
			}
		}
		sc.close();

		int[] frequencies = new int[table.size()];

		for (int i = 0; i < table.size(); i++) {

			frequencies[i] = table.get(i);
			// System.out.println(frequencies[i]);
		}
		return frequencies;
	}

	public static double[] queueToArray(Queue<Double> q) {

		int size = q.size();
		double[] b = new double[size];

		for (int k = 0; k < size; k++) {
			b[k] = q.remove();
		}
		return b;
	}

	public double[] generateSampleLambdas(String fileName) {

		int[] frequencies = readFrequencies(fileName);
		Queue<Double> queueOfLambdas = new LinkedList<Double>();
		Random random = new Random();

		double totalNumOflambdas = 0;

		for (int i = 0; i < frequencies.length; i++) {

			totalNumOflambdas += frequencies[i];
		}

		for (int i = 0; i < frequencies.length; i++) {

			int numOfRandSamplesInthisInterval = (int) ((frequencies[i] / totalNumOflambdas) * sampleSize);

			if (numOfRandSamplesInthisInterval != 0) {

				for (int j = 0; j < numOfRandSamplesInthisInterval; j++) {

					// double randomValue = rangeMin + (rangeMax - rangeMin) *
					// r.nextDouble();
					double lambda = (i * 1.0 / 10) + (0.1)
							* random.nextDouble();
					queueOfLambdas.add(lambda);
				}
			}
		}

		return queueToArray(queueOfLambdas);

	}

//	public static void main(String[] args) {
//
//		GenerateBiasedSampleOfLambdas gbsl = new GenerateBiasedSampleOfLambdas();
//
//		gbsl.sampleSize = 1000050;
//		
//		double[] lambdas = gbsl.generateSampleLambdas("/home/setas/workspace/BirthDeathModel/src/Files/lamFreqColForm");
//
//		for (int i = 0; i < lambdas.length; i++) {
//			System.out.println(lambdas[i]);
//		}
//	}
}
