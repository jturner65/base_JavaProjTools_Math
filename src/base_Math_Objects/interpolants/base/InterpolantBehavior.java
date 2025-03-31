package base_Math_Objects.interpolants.base;

import java.util.HashMap;
import java.util.Map;

/**
 * describes the behavior of interpolant when used as an animator
 * @author john
 */
public enum InterpolantBehavior {
	pingPong, 
	pingPongStop, 
	oneWayFwdLoop, 
	oneWayBkwdLoop,
	oneWayFwdStopLoop, 
	oneWayBkwdStopLoop;
	private static final String[] 
			_typeExplanation = new String[] {"Ping-pong", "Ping-pong w/pause", "1-way fwd loop", "1-way bckwd loop", "1-way fwd pause loop", "1-way bckwd pause loop"};
	private static final String[] 
			_typeName = new String[] {"Ping-pong", "Ping-pong w/pause", "1-way fwd loop", "1-way bckwd loop",  "1-way fwd pause loop", "1-way bckwd pause loop"};
	//used for file names
	private static final String[] 
			_typeBrfName = new String[] {"pingPong", "pingPongStop", "oneWayFwdLoop", "oneWayBkwdLoop", "oneWayFwdStopLoop", "oneWayBkwdStopLoop"};
	
	public static String[] getListOfTypes() {return _typeName;}
	private static Map<Integer, InterpolantBehavior> map = new HashMap<Integer, InterpolantBehavior>(); 
		static { for (InterpolantBehavior enumV : InterpolantBehavior.values()) { map.put(enumV.ordinal(), enumV);}}
	public int getVal(){return ordinal();}
	public static InterpolantBehavior getEnumByIndex(int idx){return map.get(idx);}
	public static InterpolantBehavior getEnumFromValue(int idx){return map.get(idx);}
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

}// enum InterpolantBehavior 


