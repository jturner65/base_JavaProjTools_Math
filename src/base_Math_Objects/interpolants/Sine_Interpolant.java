package base_Math_Objects.interpolants;

import base_Math_Objects.MyMathUtils;
import base_Math_Objects.interpolants.base.Base_Interpolant;

/**
 * interpolant from 1/2 period of sine
 * @author john
 *
 */
public class Sine_Interpolant extends Base_Interpolant {

    public Sine_Interpolant(float _t) {        super(_t);    }
    public Sine_Interpolant(float _t, float _stopTimer) {super(_t,_stopTimer);}

    @Override
    protected float calcInterpolant_Indiv(float _rawT) {    
        //return .5f*(1.0f + (float) Math.sin(MyMathUtils.HALF_PI *((_rawT * 2.0f) - 1.0f)));
        return .5f*(1.0f + (float) Math.sin(MyMathUtils.PI *_rawT - MyMathUtils.HALF_PI));
    }

}
