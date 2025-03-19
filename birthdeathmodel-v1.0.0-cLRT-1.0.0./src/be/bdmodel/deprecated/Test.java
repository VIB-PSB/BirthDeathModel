package be.ugent.psb.setas.bdmodel.oldclasses;

import java.text.NumberFormat;


//some testing comment
public class Test {

	public static void main(String[] args){
		String value_string = "3.826172623832312E-234";
		double value		= Double.parseDouble(value_string);
		
		double log_value	= cern.jet.math.Arithmetic.log10(value);
		
		
		NumberFormat format	= NumberFormat.getInstance();
		format.setMaximumFractionDigits(1);
		
		String new_value_string = format.format(value);
//		System.out.println(value_string);
//		System.out.println(value);		
//		System.out.println(new_value_string);
//		
//		System.out.println(log_value);
//		System.out.println(format.format(log_value));
	}
	
}