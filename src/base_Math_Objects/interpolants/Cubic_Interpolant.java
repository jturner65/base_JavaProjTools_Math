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
    protected float calcInterpolant_Indiv(float _rawT) {    return (_rawT * _rawT * (3.0f - (2.0f*_rawT)));}

}//class Hermite_Interpolant
