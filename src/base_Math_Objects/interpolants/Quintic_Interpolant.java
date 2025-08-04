package base_Math_Objects.interpolants;

import base_Math_Objects.interpolants.base.Base_Interpolant;
/**
 * 5th order interpolant with continuous accel
 * @author john
 *
 */
public class Quintic_Interpolant extends Base_Interpolant {

    public Quintic_Interpolant(float _t) {    super(_t);}
    public Quintic_Interpolant(float _t, float _stopTimer) {super(_t,_stopTimer);}

    @Override
    protected float calcInterpolant_Indiv(float _rawT) {        return (_rawT * _rawT * _rawT *(10.0f  + _rawT*(-15.0f + 6.0f*_rawT)));    }

}//class Quintic_Interpolant
