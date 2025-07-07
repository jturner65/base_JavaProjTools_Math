package base_Math_Objects.vectorObjs.floats;

import base_Math_Objects.MyMathUtils;

public class myCntlPtf extends myPointf {
    public int ID;
    public static int IDincr = 0;
    public static final float maxR = 75, 
            minR = 1,
            baseRad = 20;            //default radius for control points
    public float r, w;                //weight is calculated based on the distance to neighboring cntl myPoints when cntl myPoints are drawn
    
    public myCntlPtf(myPointf _p, float _r, float _w){ super(_p.x,_p.y, _p.z);ID=IDincr++;r=_r; w=_w; }
    public myCntlPtf(myPointf _p, float _w){this( _p, baseRad, _w);}
    public myCntlPtf(myPointf _p){this( _p, baseRad, baseRad);}
    public myCntlPtf(){this(new myPointf(), baseRad, 1);}
    //Be sure to undo IDincr increment from main constructor
    public myCntlPtf(myCntlPtf _p){this(new myPointf(_p),_p.r,_p.w);ID = _p.ID; --IDincr;}        
    public static myCntlPtf L(myCntlPtf A, float s, myCntlPtf B){    return new myCntlPtf(new myPointf(A, s, B), capInterpR(A.r, s, B.r), (1-s)*A.w + (s)*B.w);}
    public static myCntlPtf P(myCntlPtf A, myCntlPtf B){    float s = .5f;return L(A, s, B);}
    @Override
    public myPointf set(myPointf P){super.set(P); return (myPointf)this;}    
    public myCntlPtf set(myCntlPtf _p){super.set(_p.x,_p.y, _p.z); r = _p.r; w = _p.w; return this;}
    private static float capInterpR(float a, float s, float b){ float res = (1-s)*a + (s)*b; res = (res < minR ? minR : res > maxR ? maxR : res); return res;}
    
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
