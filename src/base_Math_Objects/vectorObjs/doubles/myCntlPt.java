package base_Math_Objects.vectorObjs.doubles;

public class myCntlPt extends myPoint {
	public int ID;
	public static int IDincr = 0;
	public static final float maxR = 75, 
			minR = 1,
			baseRad = 20;			//default radius for control points
	public float r, w;				//weight is calculated based on the distance to neighboring cntl myPoints when cntl myPoints are drawn
	//public static int[][] clr = new int[][]{{0,0,255,255}, {111,111,111,255}};
	//private static int[] blackClr = new int[] {0,0,0,0};
	
	public myCntlPt(myPoint _p, float _r, float _w){ super(_p.x,_p.y, _p.z);ID=IDincr++;r=_r; w=_w; }
	public myCntlPt(myPoint _p, float _w){this( _p, baseRad, _w);}
	public myCntlPt(myPoint _p){this( _p, baseRad, baseRad);}
	public myCntlPt(){this(new myPoint(),1);}
	public myCntlPt(myCntlPt _p){this(new myPoint(_p),_p.w); r = _p.r; w = _p.w;ID = _p.ID;}		
	public static myCntlPt L(myCntlPt A, float s, myCntlPt B){	return new myCntlPt(new myPoint(A, s, B), capInterpR(A.r, s, B.r), (1-s)*A.w + (s)*B.w);}//(1-s)*A.r + (s)*B.r,
	public static myCntlPt P(myCntlPt A, myCntlPt B){	float s = .5f;return L(A, s, B);}
	public myPoint set(myPoint P){super.set(P); return (myPoint)this;}
	private static float capInterpR(float a, float s, float b){ float res = (1-s)*a + (s)*b; res = (res < minR ? minR : res > maxR ? maxR : res); return res;}
	
//	public void drawMe(IRenderInterface pa, int[] clr, boolean flat){	pa.setFill(clr,clr[3]);  pa.setStroke(clr,clr[3]);		pa.show(this,2,-1,-1, flat);}		
//	public void drawRad(IRenderInterface pa, int[] clr,myVector I, myVector J){
//		pa.setFill(clr,clr[3]);  
//		pa.setStroke(clr,clr[3]);
//        pa.drawCircle(this, r, I,J,20);
//    }
//	public void drawRad(IRenderInterface pa, myVector I, myVector J){
//        pa.drawCircle(this, r, I,J,20);
//    }
//	public void drawBall(IRenderInterface pa, int[] clr,myVector I, myVector J) {
//	    float rhalf = this.r*0.5f;
//	    myPoint center1 = new myPoint(this);center1._add(myVector._mult(I,rhalf));
//	    myPoint center2 = new myPoint(this);center2._add(myVector._mult(I,-rhalf));
//	    pa.setFill(clr,clr[3]);  
//	    pa.setColorValStroke(IRenderInterface.gui_Black,255); 
//        pa.drawCircle(center1, rhalf, I,J,20);
//        pa.drawCircle(center2, rhalf, I,J,20);
//        pa.showPtAsSphere(center1, 1.0f);
//        pa.showPtAsSphere(center2, 1.0f);
//    }
//	public void drawNorm(IRenderInterface pa, int cIdx,myVector I, myVector J) {
//	    myPoint p1 = new myPoint(this);p1._add(myVector._mult(I,r));
//        myPoint p2 = new myPoint(this);p2._add(myVector._mult(I,-r));
//        pa.setColorValStroke(IRenderInterface.gui_Black,255); 
//        pa.line(p1, p2); 
//	}	
	
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

	public void calcRadFromWeight(float lenRatio, boolean inv, float wScale){r = Math.min(maxR, Math.max(minR, baseRad * (inv ? (lenRatio/w) : (wScale*w/(lenRatio*lenRatio)))));  }
	public void modRad(float modAmt){float oldR = r; r += modAmt; r = (r < minR ? minR : r > maxR ? maxR : r); w *= oldR/r; }
	public String toString(){String res = "Cntl Pt ID:"+ID+" p:"+super.toString()+" r:"+r+" w:"+w;return res;}
}//class cntlPoint