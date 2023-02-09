package base_Math_Objects.interpolants;

import base_Math_Objects.interpolants.base.Base_Interpolant;
/**
 * 5th order interpolant with continous accel
 * @author john
 *
 */
public class Quintic_Interpolant extends Base_Interpolant {

	public Quintic_Interpolant(float _t) {	super(_t);}
	public Quintic_Interpolant(float _t, float _stopTimer) {super(_t,_stopTimer);}

	@Override
	protected float calcInterpolant_Indiv(float _rawt) {		return (_rawt * _rawt * _rawt *(10.0f  + _rawt*(-15.0f + 6.0f*_rawt)));	}

}//class Quintic_Interpolant
