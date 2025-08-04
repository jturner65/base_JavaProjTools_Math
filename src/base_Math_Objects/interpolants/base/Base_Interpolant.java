package base_Math_Objects.interpolants.base;

import base_Math_Objects.interpolants.Cubic_Interpolant;
import base_Math_Objects.interpolants.Linear_Interpolant;
import base_Math_Objects.interpolants.Quintic_Interpolant;
import base_Math_Objects.interpolants.Sine_Interpolant;

/**
 * This class manages an interpolant - a value to be varied between 0 and 1
 * to be used to interpolate, extrapolate, or animate, between some set of values
 * @author john
 *
 */
public abstract class Base_Interpolant {
    
    /**
     * basic interpolant - always between 0 and 1 and linearly evolved
     */
    private float _rawT;
    /**
     * interpolant to be used by consumer - transformed 
     */
    private float t;
    
    /**
     * used to determine direction of modification for interpolant
     */
    private float modDir = 1.0f;
    
    /**
     * for interpolants that stop at extremal locations
     */
    private float stopTimeT = 0.0f;
    public final float stopTimerDur;
    public static final float _dfltStopTimerDur = 3.0f;
    private boolean isStopped;
    /**
     * used to determine desired behavior of animation when hitting extremal values 
     */
    private InterpolantBehavior animBehavior; 
    
    
    public Base_Interpolant(float _t) {
        this(_t, _dfltStopTimerDur);
    }
    
    public Base_Interpolant(float _t, float _stopTimerDur) {
        setValue(_t);
        setAnimBehavior(0);
        stopTimerDur = _stopTimerDur;    
        isStopped = false;
    }
    
    public final void setAnimBehavior(int _idx) {animBehavior = InterpolantBehavior.getEnumByIndex(_idx);}
    /**
     * Set underlying raw interpolant value - clipped to be between 0 and 1.
     * @param _t
     */
    public final void setValue(float _t) {
        _rawT = (_t<0 ? 0 : _t>1.0f ? 1.0f : _t);
        isStopped = false;    
        t = calcInterpolant(_rawT);
    }
    /**
     * return processed fade/interpolant
     * @return
     */
    public float getValue() {return t;}
    
    /**
     * Function evolves rawT, only to be used for actual animation
     * @param amt
     */
    protected final void _evolveRawT (float amt) {
        switch(animBehavior) {
            case pingPong            :
            case pingPongStop        : 
            case oneWayFwdLoop       :
            case oneWayFwdStopLoop   : {    _rawT += amt;    break;}
            case oneWayBkwdLoop      : 
            case oneWayBkwdStopLoop  : {    _rawT -= amt;    break;}            
            default                  : {    _rawT += amt;}              //unknown defaults to ping-pong
        }//switch            
    }
    
    /**
     * Reset the stop timer value and enable stop flag
     */
    private void _restartStopTime() {        stopTimeT = 0.0f;    isStopped = true; }
    
    /**
     * Evolve the stop time timer, and set whether the desired stop duration has passed
     * @param delta
     */
    private void _evolveStopTimer(float delta) {                    
        stopTimeT += delta;
        isStopped = (stopTimeT < stopTimerDur);
    }
    
    /**
     * Evolve the value of the interpolant by given amount
     * @param delta amount to evolve interpolant
     * @return
     */
    public final float evolveInterpolant(float delta) {        
        //if(isStopped) {    System.out.println("Stopped :stopTimeT : " + stopTimeT + " | delta :"+delta + " | stopTimerDur : " + stopTimerDur);}
        switch(animBehavior) {
            case pingPong         : {
                isStopped = false;
                _rawT += (modDir * delta);    
                if(_rawT >= 1.0f) {         _rawT = 1.0f; modDir = -1.0f;} 
                else if (_rawT <= 0.0f) {   _rawT = 0.0f; modDir = 1.0f;}        
                break;}
            case pingPongStop         : {
              //manage timer for duration of stop
                if(isStopped) {              _evolveStopTimer(delta);} 
                else {        
                    _rawT += (modDir * delta);    
                    if(_rawT >= 1.0f) {       _rawT = 1.0f;  modDir = -1.0f; _restartStopTime();} 
                    else if (_rawT <= 0.0f) { _rawT = 0.0f;  modDir = 1.0f;  _restartStopTime();} 
                }
                break;}
            case oneWayFwdLoop         : {
                isStopped = false;
                _rawT += delta;    
                if(_rawT > 1.0f) {     _rawT = 0.0f;}  //restart at beginning
                break;}
            case oneWayBkwdLoop     : {
                isStopped = false;
                _rawT  -= delta;    
                if (_rawT < 0.0f) {    _rawT = 1.0f;}  //restart at ending
                break;}
            case oneWayFwdStopLoop        : {
                //manage timer for duration of stop
                if(isStopped) {              _evolveStopTimer(delta);
                    if(!isStopped) {         _rawT = 0.0f;}     //stop time ended, start at beginning again
                } else {        
                    _rawT += delta;    
                    if(_rawT >= 1.0f) {      _rawT = 1.0f;      _restartStopTime();}
                }
                break;}
            case oneWayBkwdStopLoop        : {
                //manage timer for duration of stop
                if(isStopped) {              _evolveStopTimer(delta);
                    if(!isStopped) {         _rawT = 1.0f;}      //stop time ended, start at end again
                } else {        
                    _rawT  -= delta;    
                    if (_rawT <= 0.0f) {     _rawT = 0.0f;      _restartStopTime();}
                }
                break;}
            
            default         :    {//unknown defaults to ping-pong
                _rawT += (modDir * delta);    
                if(_rawT > 1.0f) {_rawT = 1.0f;modDir = -1.0f;} else if (_rawT < 0.0f) {    _rawT = 0.0f;    modDir = 1.0f;}                    
            }
        }//switch
//        _rawT += (modDir * delta);            
//        if(_rawT > 1.0f) {_rawT = 1.0f;modDir = -1.0f;} else if (_rawT < 0.0f) {    _rawT = 0.0f;    modDir = 1.0f;}
        t = calcInterpolant(_rawT);
        return t;
    }
    /**
     * build actual interpolant fade 
     * @return
     */
    public final float calcInterpolant(float _rawT) {
        if(_rawT<=0) {return 0.0f;}
        if(_rawT>=1.0f) {return 1.0f;}
        return calcInterpolant_Indiv(_rawT);
    }
    protected abstract float calcInterpolant_Indiv(float _rawT);
    
    /**
     * build an interpolant of passed type, for type defined in Base_Interpolant
     * @param animType
     * @return
     */
    public static Base_Interpolant buildInterpolant(InterpolantTypes animType, float _initT) {return buildInterpolant(animType, _initT, _dfltStopTimerDur);}
    public static Base_Interpolant buildInterpolant(InterpolantTypes animType, float _initT, float _stopTime) {
        switch(animType) {
            case linear          : {        return new Linear_Interpolant(_initT,_stopTime);        }
            case smoothVelocity  : {        return new Cubic_Interpolant(_initT,_stopTime);        }
            case smoothAccel     : {        return new Quintic_Interpolant(_initT,_stopTime);        }
            case sine            : {        return new Sine_Interpolant(_initT,_stopTime);        }
            default : {
                System.out.println("Base_Interpolant :: buildInterpolant :: Unknown interpolant type : " + animType.toString() + ".  Aborting.");
                return null;
            }
        }
    }//buildInterpolant
    

}//class Base_Interpolant

