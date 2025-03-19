package be.ugent.psb.setas.bdmodel.model;

public class DoubleTesting {


		static final double EPSILON = 1e-300;
		
		public static boolean isEqual(double a, double b, double PRECISION){
			return (Math.abs(a-b) < PRECISION);
		}
		public static void main(String[] args){
			
//			System.out.println("MIN_EXPONENT: " + Double.MIN_EXPONENT);
//			System.out.println("MIN_NORMAL: " + Double.MIN_NORMAL);
//			System.out.println("MIN_VALUE: " + Double.MIN_VALUE);
//			
//			System.out.println();
//			
//			System.out.println("MAX_EXPONENT: " + Double.MAX_EXPONENT);
//			System.out.println("MAX_VALUE: " + Double.MAX_VALUE);
			
			
			
			double a = 1e-250;
			double b = 2e-250;
			
			
			if(a < b){
				System.out.println("A < B");
			}
			if (isEqual(a,b, EPSILON)){
				System.out.println("A == B");
			} else {
				System.out.println("A =|= B");
			}
			
			
			
		}
		
		
	}

