package base_Math_Objects.interpolants;

import base_Math_Objects.interpolants.base.Base_Interpolant;

/**
 * hermitian blending func with continuous vel
 * @author john
 *
 */
public class Cubic_Interpolant extends Base_Interpolant {

	public Cubic_Interpolant(float _t) {super(_t);}
	public Cubic_Interpolant(float _t, float _stopTimer) {super(_t,_stopTimer);}

	@Override
	protected float calcInterpolant_Indiv(float _rawt) {	return (_rawt * _rawt * (3.0f - (2.0f*_rawt)));}

}//class Hermite_Interpolant
