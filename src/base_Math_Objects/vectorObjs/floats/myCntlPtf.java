package base_Math_Objects.vectorObjs.floats;

import base_Math_Objects.MyMathUtils;

/**
 * Specialization of a floating-point based myPointf intended to be used as a control point for trajectories.
 */
public class myCntlPtf extends myPointf {
    /**
     * Control point ID
     */
    public int ID;
    private static int IDincr = 0;
    /**
     * Radius values for display of control points
     */
    public static final float 
            maxR = 75, 
            minR = 1,
            baseRad = 20;            //default radius for control points
    public float r, w;                //weight is calculated based on the distance to neighboring cntl myPoints when cntl myPoints are drawn
    
    public myCntlPtf(myPointf _p, float _r, float _w){ super(_p.x, _p.y, _p.z);ID=IDincr++;r=_r; w=_w; }
    public myCntlPtf(myPointf _p, float _w){this( _p, baseRad, _w);}
    public myCntlPtf(myPointf _p){this( _p, baseRad, baseRad);}
    public myCntlPtf(){this(new myPointf(), baseRad, 1);}
    /**
     * Copy ctor
     * @param _p
     */
    public myCntlPtf(myCntlPtf _p){ super(_p);  ID = _p.ID;  r = _p.r; w = _p.w;}
    /**
     * Interpolating constructor
     * @param A
     * @param s
     * @param B
     */
    public myCntlPtf(myCntlPtf A, float s, myCntlPtf B) {
        super(A, s, B);
        r = _cappedLinInterp(A.r, s, B.r, minR, maxR);
        w = _linInterp(A.w, s, B.w);
    }    
    /**
     * Interpolating constructor - equidistant between both passed points 
     * @param A
     * @param B
     */
    public myCntlPtf(myCntlPtf A, myCntlPtf B) {     this(A, .5f, B);  }
    @Override
    public myPointf set(myPointf P){super.set(P); return (myPointf)this;}    
    public myCntlPtf set(myCntlPtf _p){super.set(_p.x,_p.y, _p.z); r = _p.r; w = _p.w; return this;}
    /**
     * calc the rotation of this point by angle a around G on plane described by I, in direction inferred by J(tangent)
     * @param a
     * @param I
     * @param J
     * @param G
     * @return
     */
    @Override
    public final myCntlPtf rotMeAroundPt(float a, myVectorf I, myVectorf J, myPointf G) {
        return new myCntlPtf(super.rotMeAroundPt(a, I, J, G), r, w); 
    }; 
    
    
    /**
     * returns rotated version of this point by angle(CP,CR) parallel to plane (C,P,R)
     * @param C
     * @param P
     * @param R
     * @return this point rotated
     */
    @Override
    public final myCntlPtf rotMeAroundPt(myPointf C, myPointf P, myPointf R) { // returns rotated version of Q by angle(CP,CR) parallel to plane (C,P,R) 
        return new myCntlPtf(super.rotMeAroundPt(C,P,R), r, w); 
    }    
    
    public void calcRadFromWeight(float lenRatio, boolean inv, float wScale){r = MyMathUtils.min(maxR, MyMathUtils.max(minR, baseRad * (inv ? (lenRatio/w) : (wScale*w/(lenRatio*lenRatio)))));  }
    public void modRad(float modAmt){float oldR = r; r += modAmt; r = (r < minR ? minR : r > maxR ? maxR : r); w *= oldR/r; }
    @Override
    public String toString(){String res = "Cntl Pt_f ID:"+ID+" p:"+super.toString()+" r:"+r+" w:"+w;return res;}
}
