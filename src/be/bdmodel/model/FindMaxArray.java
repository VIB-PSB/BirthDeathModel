package be.ugent.psb.setas.bdmodel.model;

public class FindMaxArray {

	public static double[] FindMax(double[] a) {
		
		double maxVal = a[1];
		double maxIndice = 1;
		
		double[] A = new double[2];
		for (int i = 2; i < a.length; i++) {
			if (a[i] > maxVal) {
				maxVal = a[i];
				maxIndice = i;
			}
		}
		A[0] = maxVal;
		A[1] = maxIndice;
		return A;
	}

	public static double findMaxValue(double[] a) {
		//because in this be.ugent.psb.setas.bdmodel.model we do not deal with rootSize =0, and the likelihoods are in log scale <0
		double maxVal = a[1];
		for (int i = 2; i < a.length; i++) {
			if (a[i] > maxVal) {
				maxVal = a[i];
			}
		}
		return maxVal;
	}

	public static int findMaxValueIntArray(int[] a) {
		int maxVal = a[1];
		for (int i = 2; i < a.length; i++) {
			if (a[i] > maxVal) {
				maxVal = a[i];
			}
		}
		return maxVal;
	}

	public static int findIndexOfMaxValue(double[] a) {
		int index = 1;
		for (int i = 2; i < a.length; i++) {
			if (a[i] > a[index]) {
				index = i;
			}
		}
		return index;
	}

	
	public static int[] findIndexOfMaxValueInMatrix(double[][] a, int rows, int cols) {
		
		int [] index =  new int [2];
		index[0] =1;
		index[1]=1;
		
		for (int i = 1; i < rows; i++) {
			for(int j=0; j< cols;j++){
			if (a[i][j] > a[index[0]][index[1]]) {
				index[0] = i;
				index[1]=j;
			}
		}
		}
		return index;
	}
	
	public static int findIndexOfMaxValueIntArray(int[] a) {
		int index = 1;
		for (int i = 1; i < a.length; i++) {
			if (a[1] > a[index]) {
				index = i;
			}
		}
		return index;
	}

	 public static void main(String[] args) {
	 // TODO Auto-generated method stub
	 FindMaxArray fma = new FindMaxArray();
	 double []a = new double [5];
	 a[0]= 1.1;
	 a[1] = 2.1;
	 a[2] =3.1;
	 
	 
	 a[3] =4.1;
	 a[4]= 5.1;
	 double [] b = new double [2];
	 b = fma.FindMax(a);
	 // double b = fma.ValueFindMax(a);
	 // double d = fma.IndiceFindMax(a);
	 for (int i=0; i<2;i++){
	 System.out.print(b[i]+"\t");}
	 System.out.print("\n");
	 // System.out.println(b);
	 // System.out.println(d);
	 }

}
