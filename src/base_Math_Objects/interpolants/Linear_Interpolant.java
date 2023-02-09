package base_Math_Objects.interpolants;

import base_Math_Objects.interpolants.base.Base_Interpolant;

/**
 * just linearly vary t from 0->1
 * @author john
 *
 */
public class Linear_Interpolant extends Base_Interpolant {

	public Linear_Interpolant(float _t) {	super(_t);}
	public Linear_Interpolant(float _t, float _stopTimer) {super(_t,_stopTimer);}

	@Override
	protected float calcInterpolant_Indiv(float _rawT) {return _rawT;}

}//class Linear_Interpolant
