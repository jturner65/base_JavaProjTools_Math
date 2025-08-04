package base_Math_Objects.vectorObjs.doubles;

import base_Math_Objects.MyMathUtils;

/**
 * Specialization of a double based myPoint intended to be used as a control point for trajectories.
 */
public class myCntlPt extends myPoint {
    /**
     * Control point ID
     */
    public int ID;
    public static int IDincr = 0;
    /**
     * Radius values for display of control points
     */
    public static final double maxR = 75, 
            minR = 1,
            baseRad = 20;            //default radius for control points
    public double r, w;                //weight is calculated based on the distance to neighboring cntl myPoints when cntl myPoints are drawn
    
    public myCntlPt(myPoint _p, double _r, double _w){ super(_p.x,_p.y, _p.z);ID=IDincr++;r=_r; w=_w; }
    public myCntlPt(myPoint _p, double _w){this( _p, baseRad, _w);}
    public myCntlPt(myPoint _p){this( _p, baseRad, baseRad);}
    public myCntlPt(){this(new myPoint(), baseRad, 1);}
    /**
     * Copy ctor
     * @param _p
     */
    public myCntlPt(myCntlPt _p){ super(_p);  ID = _p.ID;  r = _p.r; w = _p.w;}
    /**
     * Interpolating constructor
     * @param A
     * @param s
     * @param B
     */
    public myCntlPt(myCntlPt A, double s, myCntlPt B) {
        super(A, s, B);
        r = _cappedLinInterp(A.r, s, B.r, minR, maxR);
        w = _linInterp(A.w, s, B.w);
    }    
    /**
     * Interpolating constructor - equidistant between both passed points 
     * @param A
     * @param B
     */
    public myCntlPt(myCntlPt A, myCntlPt B) {     this(A, .5, B);  }
    @Override
    public myPoint set(myPoint P){super.set(P); return (myPoint)this;}
    public myCntlPt set(myCntlPt _p){super.set(_p.x,_p.y, _p.z); r = _p.r; w = _p.w; return this;}
    
    /**
     * calc the rotation of this point by angle a around G on plane described by I, in direction inferred by J(tangent)
     * @param a
     * @param I
     * @param J
     * @param G
     * @return
     */
    @Override
    public final myCntlPt rotMeAroundPt(double a, myVector I, myVector J, myPoint G) {
        return new myCntlPt(super.rotMeAroundPt(a, I, J, G), r, w); 
    }; 
        
    /**
     * returns rotated version of this point by angle(CP,CR) parallel to plane (C,P,R)
     * @param C
     * @param P
     * @param R
     * @return this point rotated
     */
    @Override
    public final myPoint rotMeAroundPt(myPoint C, myPoint P, myPoint R) { // returns rotated version of Q by angle(CP,CR) parallel to plane (C,P,R) 
        return new myCntlPt(super.rotMeAroundPt(C,P,R), r, w); 
    }     

    public void calcRadFromWeight(double lenRatio, boolean inv, double wScale){r = MyMathUtils.min(maxR, MyMathUtils.max(minR, baseRad * (inv ? (lenRatio/w) : (wScale*w/(lenRatio*lenRatio)))));  }
    public void modRad(double modAmt){double oldR = r; r += modAmt; r = (r < minR ? minR : r > maxR ? maxR : r); w *= oldR/r; }
    public String toString(){String res = "Cntl Pt ID:"+ID+" p:"+super.toString()+" r:"+r+" w:"+w;return res;}
}//class cntlPoint