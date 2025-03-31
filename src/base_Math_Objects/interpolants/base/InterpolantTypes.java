package base_Math_Objects.interpolants.base;

import java.util.HashMap;
import java.util.Map;


/**
 * describes the type of interpolant used
 * @author john
 *
 */
public enum InterpolantTypes {	
	linear, smoothVelocity, smoothAccel, sine;
	private static final String[] 
			_typeExplanation = new String[] {"Linear","Cubic (Continuous Velocity)","Quintic (Continuous Accel)","Sine"};
	private static final String[] 
			_typeName = new String[] {"Linear","Cubic","Quintic","Sine"};
	//used for file names
	private static final String[] 
			_typeBrfName = new String[] {"linear","cubic","quintic","sine"};
	
	public static String[] getListOfTypes() {return _typeName;}
	private static Map<Integer, InterpolantTypes> map = new HashMap<Integer, InterpolantTypes>(); 
		static { for (InterpolantTypes enumV : InterpolantTypes.values()) { map.put(enumV.ordinal(), enumV);}}
	public int getVal(){return ordinal();}
	public static InterpolantTypes getVal(int idx){return map.get(idx);}
	public static int getNumVals(){return map.size();}						//get # of values in enum
	public String getName() {return _typeName[ordinal()];}
	public String getBrfName() {return _typeBrfName[ordinal()];}
	public String getExplanation() {return _typeExplanation[ordinal()];}
	public static String getNameByVal(int _val) {return _typeName[_val];}
	public static String getBrfNameByVal(int _val) {return _typeBrfName[_val];}
	public static String getExplanationByVal(int _val) {return _typeExplanation[_val];}
	@Override
    public String toString() { return ""+this.name()+":"+_typeExplanation[ordinal()]; }	
    public String toStrBrf() { return ""+_typeExplanation[ordinal()]; }	

}


