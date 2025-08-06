package base_Math_Objects.vectorObjs.doubles;

import base_Math_Objects.vectorObjs.floats.myPointf;

public class myPoint {
    /**
     * this point object's x coord
     */
    public double x;
    /**
     * this point object's y coord
     */
    public double y;
    /**
     * this point object's z coord
     */
    public double z;
    /**
     * Static final member : origin point in 3D
     */
    public static final myPoint ZEROPT = new myPoint(0,0,0);

    /**
     * build a point with given coordinates
     * @param _x : x coord
     * @param _y : y coord
     * @param _z : z coord
     */
    public myPoint(double _x, double _y, double _z){this.x = _x; this.y = _y; this.z = _z;}         //constructor 3 args
    
    /**
     * Build a point using the first 3 values in the passed array. Will fail if less than 3 values.
     * @param vals must be 3+ values in length. Any values past the first 3 will be ignored.
     */
    public myPoint(float[] vals){this.x = vals[0]; this.y = vals[1]; this.z = vals[2];} 
    
    /**
     * Build a point using the first 3 values in the passed array. Will fail if less than 3 values.
     * @param vals must be 3+ values in length. Any values past the first 3 will be ignored.
     */
    public myPoint(double[] vals){this.x = vals[0]; this.y = vals[1]; this.z = vals[2];} 
    /**
     * copy constructor
     * @param p : point object to copy
     */
    public myPoint(myPoint p){ this(p.x, p.y, p.z); }                                                                     
    /**
     * build point as displacement from point A by vector B
     * @param A : starting point
     * @param B : displacement vector
     */
    public myPoint(myPoint A, myVector B) {this(A.x+B.x,A.y+B.y,A.z+B.z); };
    /**
     * Interpolate between A and B by s -> (0->1)
     * @param A : first point to interpolate from
     * @param s : value [0,1] to determine linear interpolation
     * @param B : second point to interpolate from
     */
    public myPoint(myPoint A, double s, myPoint B) {   this(_linInterp(A.x, s, B.x), _linInterp(A.y, s, B.y), _linInterp(A.z, s, B.z)); };        //builds a point somewhere in between a and b
    /**
     * empty constructor
     */
    public myPoint(){ this(0,0,0);}                                                                                                                               //constructor 0 args
    /**
     * clear this point's set values
     */
    public void clear() {this.x = 0; this.y = 0; this.z = 0;}    
    /**
     * Set this object's coordinate values
     * @param _x, _y, _z : new x,y,z coords of this object
     */
    public void set(double _x, double _y, double _z){ this.x = _x;  this.y = _y;  this.z = _z; }                                               //set 3 args 
    /**
     * Set this point's values as a copy of the passed point
     * @param p : the point to copy
     * @return this point
     */
    public myPoint set(myPoint p){ this.x = p.x; this.y = p.y; this.z = p.z; return this;}                                                                   //set 1 args
    /**
     * build and return the average (midpoint) of this point and the passed point
     * @param q : the point to find the midpoint with
     * @return the midpoint between this and q
     */
    public myPoint _avgWithMe(myPoint q) {return new myPoint((this.x+q.x)/2.0,(this.y+q.y)/2.0,(this.z+q.z)/2.0);} 
    /**
     * build and return the average (midpoint) of this point and the passed points
     * @param q,r : the points to find the midpoint with
     * @return the midpoint between this and q and r
     */    
    public myPoint _avgWithMe(myPoint q, myPoint r) {return new myPoint((this.x+q.x+r.x)/3.0f,(this.y+q.y+r.y)/3.0f,(this.z+q.z+r.z)/3.0f);} 
    /**
     * build and return the average (midpoint) of this point and the passed points
     * @param q,r : the points to find the midpoint with
     * @return the midpoint between this and q and r
     */    
    public myPoint _avgWithMe(myPoint q, myPoint r, myPoint s) {return new myPoint((this.x+q.x+r.x+s.x)/4.0f,(this.y+q.y+r.y+s.y)/4.0f,(this.z+q.z+r.z+s.z)/4.0f);}
    /**
     * Static method : build the midpoint between the two passed points
     * @param p,q : 2 points to find the midpoint between
     * @return midpoint between p and q
     */
    public static myPoint _average(myPoint p, myPoint q) {return new myPoint((p.x+q.x)/2.0,(p.y+q.y)/2.0,(p.z+q.z)/2.0);}     
    /**
     * Static method : build the midpoint between the three passed points
     * @param p,q,r : 3 points to find the midpoint between
     * @return midpoint between  p,q,r
     */
    public static myPoint _average(myPoint p, myPoint q, myPoint r) {return new myPoint((p.x+q.x+r.x)/3.0f,(p.y+q.y+r.y)/3.0f,(p.z+q.z+r.z)/3.0f);} 
    /**
     * Static method : build the midpoint between the two passed points
     * @param p,q,r,s : 4 points to find the midpoint between
     * @return midpoint between p,q,r,s
     */
    public static myPoint _average(myPoint p, myPoint q, myPoint r, myPoint s) {return new myPoint((p.x+q.x+r.x+s.x)/4.0f,(p.y+q.y+r.y+s.y)/4.0f,(p.z+q.z+r.z+s.z)/4.0f);} 
    /**
     * Static method : build the center of verticies for the array of passed points (COV)
     * @param pts : array of points
     * @return COV for pts array
     */
    public static myPoint _average(myPoint[] pts) {myPoint res = new myPoint();    for(myPoint p : pts) {    res._add(p);}return myPoint._div(res, 1.0f*pts.length);}     /**
     * multiply this point by n, modifying this point to be result
     * @param n : value to scale this vector by
     * @return this point, after scaling
     */
    public myPoint _mult(double n){ this.x *= n; this.y *= n; this.z *= n; return this; }                                                     //_mult 3 args  
    /**
     * Static method : return result of multiplying a point by a scalar
     * @param p : point to be scaled
     * @param n : scale value (to multiply)
     * @return point result of element-wise multiplication of p by n
     */
    public static myPoint _mult(myPoint p, double n){ return new myPoint(p.x * n, p.y * n, p.z * n);}                          //1 pt, 1 double
    /**
     * Static method : element-wise multiplication of two points, returning result
     * @param p,q : two points to multiply element-wise with each other
     * @return : element-wise product of p and q : (p.x * q.x) i + (p.y*q.y) j + (p.z*q.z) k
     */
    public static myPoint _mult(myPoint p, myPoint q){ return new myPoint(p.x *q.x, p.y * q.y, p.z * q.z);}           //return elementwise product
    /**
     * Static method : element-wise multiplication of two points, returning result
     * @param p,q : two points to element-wise multiply
     * @param r : destination for result of element-wise multiplication
     */
    public static void _mult(myPoint p, myPoint q, myPoint r){ myPoint result = new myPoint(p.x *q.x, p.y * q.y, p.z * q.z); r.set(result);}           //2 pt src, 1 pt dest      
    /**
     * divide this point by q, making this point equal to result.  No value checking is performed
     * @param q scalar value to divide this point by
     */
    public void _div(double q){this.x /= q; this.y /= q; this.z /= q; }  
    /**
     * Static Method : divide passed point p by n, making this point equal to result.  No value checking is performed
     * @param p : point to scale (divide)
     * @param n : value to divide p by
     * @return point result of per-element division of p by n
     */
    public static myPoint _div(myPoint p, double n){ if(n==0) return p; return new myPoint(p.x / n, p.y / n, p.z / n);}                          //1 pt, 1 double
    /**
     * add passed values to this point, making this point equal to result. 
     * @param _x : x coord to add to this.x
     * @param _y : y coord to add to this.y
     * @param _z : z coord to add to this.z
     */
    public void _add(double _x, double _y, double _z){ this.x += _x; this.y += _y; this.z += _z;   }  
    /**
     * element-wise add passed point to this point, making this point equal to result. 
     * @param v point to add to this
     */
    public void _add(myPoint v){ this.x += v.x; this.y += v.y; this.z += v.z;   }                                                 //_add 1 arg  
    /**
     * Static Method : Add vector I to point O, returning result
     * @param O : origin point
     * @param I : displacement vector
     * @return : O + I
     */
    public static myPoint _add(myPoint O, myVector I){                                                            return new myPoint(O.x+I.x,O.y+I.y,O.z+I.z);}  
    /**
     * Static Method : Add vector I (scaled by a) to point O, returning result
     * @param O : origin point
     * @param a : scaling of I
     * @param I : displacement vector
     * @return : O + aI
     */
    public static myPoint _add(myPoint O, double a, myVector I){                                                return new myPoint(O.x+a*I.x,O.y+a*I.y,O.z+a*I.z);}                                        //2 vec
    /**
     * Static Method : Add vector I (scaled by a) and vector J (scaled by b) to point O, returning result
     * @param O : origin point
     * @param a : scaling of I
     * @param I : displacement vector
     * @param b : scaling of J
     * @param J : displacement vector
     * @return : O + aI + bJ
     */
    public static myPoint _add(myPoint O, double a, myVector I, double b, myVector J) {                            return new myPoint(O.x+a*I.x+b*J.x,O.y+a*I.y+b*J.y,O.z+a*I.z+b*J.z);}                      // O+xI+yJ
    /**
     * Static Method : Add vector I (scaled by a), vector J (scaled by b) and vector K (scaled by c) to point O, returning result
     * @param O : origin point
     * @param a : scaling of I
     * @param I : displacement vector
     * @param b : scaling of J
     * @param J : displacement vector
     * @param c : scaling of K
     * @param K : displacement vector
     * @return : O + aI + bJ + cK
     */
    public static myPoint _add(myPoint O, double a, myVector I, double b, myVector J, double c, myVector K) {    return new myPoint(O.x+a*I.x+b*J.x+c*K.x,O.y+a*I.y+b*J.y+c*K.y,O.z+a*I.z+b*J.z+c*K.z);} // O+xI+yJ+kZ
    /**
     * Static Method : add two points and return result
     * @param p,q : points to add
     * @return resulting point of element-wise addition of p + q
     */
    public static myPoint _add(myPoint p, myPoint q){return new myPoint(p.x + q.x, p.y + q.y, p.z + q.z); }
    /**
     * Static Method : add two points, putting result in 3rd argument
     * @param p,q : points to add
     * @param r : resulting point of element-wise addition of p + q
     */
    public static void _add(myPoint p, myPoint q, myPoint r){ myPoint result = new myPoint(p.x + q.x, p.y + q.y, p.z + q.z); r.set(result);}           //2 pt src, 1 pt dest  
    /**
     * Static Method : add an array of points, returning result
     * @param pAra
     * @return
     */
    public static myPoint _add(myPoint[] pAra) {
        // aggregate to minimize magnitude calcs
        float _x = 0, _y = 0, _z = 0;
        for(int i=0;i<pAra.length;++i) {
            _x += pAra[i].x;
            _y += pAra[i].y;
            _z += pAra[i].z;
        }
        return new myPoint(_x,_y,_z);
    }//_add

    /**
     * subtract passed values to this point, making this point equal to result. 
     * @param _x : x coord to subtract from this.x
     * @param _y : y coord to subtract from this.y
     * @param _z : z coord to subtract from this.z
     */
    public void _sub(double _x, double _y, double _z){ this.x -= _x; this.y -= _y; this.z -= _z;  }                                                                   //_sub 3 args
    /**
     * element-wise subtract passed point from this point, making this point equal to result. 
     * @param v point to subtract this
     */    
    public void _sub(myPoint v){ this.x -= v.x; this.y -= v.y; this.z -= v.z;  }                                                                           //_sub 1 arg 
    /**
     * Static Method : subtract two points and return result
     * @param p,q : points to subtract
     * @return resulting point of element-wise subtraction of p - q
     */
    public static myPoint _sub(myPoint p, myPoint q){ return new myPoint(p.x - q.x, p.y - q.y, p.z - q.z);}
    /**
     * Static Method : subtract two points, putting result in 3rd argument
     * @param p,q : points to subtract
     * @param r : resulting point of element-wise subtraction of p - q
     */
    public static void _sub(myPoint p, myPoint q, myPoint r){ myPoint result = new myPoint(p.x - q.x, p.y - q.y, p.z - q.z); r.set(result);}       //2 pt src, 1 pt dest      
    
    /**
     * create a new copy point of this point
     * @return
     */
    public myPoint cloneMe(){return new myPoint(this.x, this.y, this.z); }  
    
    /**
     * calculate L1 (Manhattan) distance between this point and passed point.  Manhattan distance is
     * Math.abs((this.x - q.x)) + Math.abs((this.y - q.y)) + Math.abs((this.z - q.z))
     * @param q : point to find distance from
     * @return : L1 distance from this point to q 
     */    
    public double _L1Dist(myPoint q){return Math.abs((this.x - q.x)) + Math.abs((this.y - q.y)) + Math.abs((this.z - q.z)); }
    /**
     * Static Method : calculate L1 (Manhattan) distance between passed points.  Manhattan distance is
     * Math.abs((r.x - q.x)) + Math.abs((r.y - q.y)) + Math.abs((r.z - q.z))
     * @param q,r : points to find distance between
     * @return : L1 distance from q to r
     */    
    public static double _L1Dist(myPoint q, myPoint r){ return q._L1Dist(r);}
    
    /**
     * find squared L2 (Euclidean) distance from this point to q.  Squared L2 distance is
     * ((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z)) 
     * @param q point to find distance from
     * @return squared L2 Distance from this point to q
     */
    public double _SqrDist(myPoint q){ double dx=(this.x-q.x), dy=(this.y-q.y), dz=(this.z-q.z);return ((dx*dx) + (dy*dy) + (dz*dz)); }
    /**
     * find squared L2 (Euclidean) distance from this point to q (point represented as float).  Squared L2 distance is
     * ((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z)) 
     * @param q point to find distance from
     * @return squared L2 Distance from this point to q
     */
    public double _SqrDist(myPointf q){ double dx=(this.x-q.x), dy=(this.y-q.y), dz=(this.z-q.z);return ((dx*dx) + (dy*dy) + (dz*dz)); }
    /**
     * Static Method : find squared L2 (Euclidean) distance from point q to point r.  Squared L2 distance is
     * ((r.x - q.x)*(r.x - q.x)) + ((r.y - q.y)*(r.y - q.y)) + ((r.z - q.z)*(r.z - q.z)) 
     * @param q,r : points to find distance between
     * @return squared L2 Distance from q to r
     */
    public static double _SqrDist(myPoint q, myPoint r){ double dx=(r.x-q.x), dy=(r.y-q.y), dz=(r.z-q.z);return ((dx*dx) + (dy*dy) + (dz*dz)); }    
    /**
     * find L2 (Euclidean) distance from this point to q.  L2 (Euclidean) distance is
     * sqrt(((r.x - q.x)*(r.x - q.x)) + ((r.y - q.y)*(r.y - q.y)) + ((r.z - q.z)*(r.z - q.z))) 
     * @param q point to find distance to
     * @return L2 Distance from this point to q
     */
    public double _dist(myPoint q){ double dx=(this.x-q.x), dy=(this.y-q.y), dz=(this.z-q.z);return Math.sqrt(((dx*dx) + (dy*dy) + (dz*dz)));}
    /**
     * Static Method : find L2 (Euclidean) distance from point q to point r.  Squared L2 distance is
     * sqrt(((r.x - q.x)*(r.x - q.x)) + ((r.y - q.y)*(r.y - q.y)) + ((r.z - q.z)*(r.z - q.z))) 
     * @param q,r : points to find distance between
     * @return L2 Distance from q to r
     */
    public static double _dist(myPoint q, myPoint r){double dx=(r.x-q.x), dy=(r.y-q.y), dz=(r.z-q.z);return Math.sqrt(((dx*dx) + (dy*dy) + (dz*dz)));}
    /**
     * find L2 (Euclidean) distance from this point to passed coordinates.  Squared L2 distance is
     * sqrt(((this.x - qx)*(this.x - qx)) + ((this.y - qy)*(this.y - qy)) + ((this.z - qz)*(this.z - qz)))
     * @param qx,qy,qz : coordinates to find distance to
     * @return L2 Distance from this to [qx,qy,qz]
     */
    public double _dist(double qx, double qy, double qz){ double dx=(this.x-qx), dy=(this.y-qy), dz=(this.z-qz);return Math.sqrt(((dx*dx) + (dy*dy) + (dz*dz)));}
    /**
     * Static Method : find L2 (Euclidean) distance from point q to passed coordinates.  Squared L2 distance is
     * sqrt(((r.x - qx)*(r.x - qx)) + ((r.y - qy)*(r.y - qy)) + ((r.z - qz)*(r.z - qz)))
     * @param r : point to find distance from
     * @param qx,qy,qz : coordinates to find distance to
     * @return L2 Distance from r to [qx,qy,qz]
     */
    public static double _dist(myPoint r, double qx, double qy, double qz){ double dx=(r.x-qx), dy=(r.y-qy), dz=(r.z-qz);return Math.sqrt((dx*dx) + (dy*dy) + (dz*dz)); }  
    /**
     * return the values of this point as an array of doubles
     * @return array of doubles {x,y,z}
     */
    public double[] asArray(){return new double[]{x,y,z};}
    /**
     * return the values of this point as a float
     * @return array of float {x,y,z}
     */    
    public float[] asFltArray(){return new float[]{(float)x,(float)y,(float)z};}
    /**
     * return the values of this point as a homogenous point array
     * @return array of doubles {x,y,z, 1}
     */
    public double[] asHAraPt(){return new double[]{this.x, this.y, this.z,1};}
    /**
     * return the values of this point as a homogenous vector array
     * @return array of doubles {x,y,z, 0}
     */
    public double[] asHAraVec(){return new double[]{this.x, this.y, this.z,0};}
    
    /**
     * this will calculate normalized barycentric coordinates (u, v, w) for pt w/respect to first 3 control points
     * @param cntlPts control points of poly
     * @return
     */
    public final double[] calcNormBaryCoords(myPoint[] cntlPts) {        
        //pt = u * cntlPts[0] + v * cntlPts[1] + w * cntlPts[2]
        myVector AB = new myVector(cntlPts[0],cntlPts[1]),
                AC = new myVector(cntlPts[0],cntlPts[2]),
                AP = new myVector(cntlPts[0],this);
        double d00 = AB.sqMagn, d01 = AB._dot(AC), 
            d11 = AC.sqMagn, d20 = AP._dot(AB), d21 = AP._dot(AC);
        double 
            denom = d00 * d11 - d01 * d01,
            v = (d11 * d20 - d01 * d21) / denom,
            w =  (d00 * d21 - d01 * d20) / denom,
            u = 1.0f - v- w;        
        return new double[] {u,v,w};
    }

    /**
     * Returns array of distances from first point in given point trajectory array to each subsequent point. 
     * @param pts array of points representing a trajectory. 
     * @param wrap if true then last entry is length of entire closed loop back to first point, otherwise length of entire curve to final point.
     * @return array of distances from first point to each of subsequent points
     */
    public static final double[] _findAllTrajPtDists(myPoint[] pts, boolean wrap){
        double[] res = new double[pts.length+1];
        res[0]=0;
        for(int i=1; i<pts.length; ++i){res[i] = res[i-1] + _dist(pts[i-1],pts[i]);}
        if(wrap){res[pts.length] = res[pts.length-1] + _dist(pts[pts.length-1],pts[0]); } 
        else {res[pts.length] = res[pts.length-1];}        
        return res;
     }
    /**
     * Returns length of curve described by array of points
     * @param pts array of points describing the curve we want the length of.
     * @param closed whether the array is a closed loop or not
     * @return the length of either the closed loop or open trajectory described by the given array of points.
     */
    public static final double _getLengthOfTraj(myPoint[] pts, boolean closed){
        double res = 0;
        for(int i=0;i<pts.length-1;++i) {res += _dist(pts[i],pts[i+1]);}
        if(closed) {res += _dist(pts[pts.length-1], pts[0]);}
        return res;
    }
    
    /**
     * given control points and passed normalized barycentric coordinates, calculate resultant point
     * @param cntlPts
     * @param pointNBC
     * @return
     */    
    public static final myPoint _calcPointFromNormBaryCoords(myPoint[] cntlPts, double[] pointNBC) {
        myPoint res = new myPoint();
        for(int i=0;i<pointNBC.length;++i) {    res._add(myPoint._mult(cntlPts[i], pointNBC[i]));}    
        return res;
    }
    
    /**
     * calc the rotation of this point by angle a around G on plane described by I, in direction inferred by J(tangent)
     * @param a
     * @param I
     * @param J
     * @param G
     * @return
     */
    public myPoint rotMeAroundPt(double a, myVector I, myVector J, myPoint G) {
        double x= myVector._dot(new myVector(G,this),myVector._unit(I)), y=myVector._dot(new myVector(G,this),myVector._unit(J)); 
        double c=Math.cos(a), s=Math.sin(a); 
        // subtract G == (x, y) from result for translation to origin 
        double iXVal = x*c-x-y*s, jYVal= x*s+y*c-y;            
        return myPoint._add(this,iXVal,I,jYVal,J); 
    }
    
    /**
     * returns rotated version of this point by angle(CP,CR) parallel to plane (C,P,R)
     * @param C
     * @param P
     * @param R
     * @return this point rotated
     */
    public myPoint rotMeAroundPt(myPoint C, myPoint P, myPoint R) { // returns rotated version of this point by angle(CP,CR) parallel to plane (C,P,R)
        myVector I0=myVector._unit(C,P), I1=myVector._unit(C,R), V=new myVector(C,this); 
        double c=myPoint._dist(I0,I1), s=Math.sqrt(1.-(c*c)); 
        if(Math.abs(s)<0.00001) return this;        
        myVector J0=myVector._add(myVector._mult(I1,1./s),myVector._mult(I0,-c/s));  
        myVector J1=myVector._add(myVector._mult(I0,-s),myVector._mult(J0,c));  
        double x=V._dot(I0), y=V._dot(J0);  
        return myPoint._add(this,x,myVector._sub(I1,I0),y,myVector._sub(J1,J0)); 
    }     
    
    /**
     * Returns whether the passed point is within eps of this point
     * @param p
     * @param eps
     * @return
     */    
    public boolean clickIn(myPoint p, double eps) { return(_dist(p) < eps);}
    
    /**
     * returns if this myPoint is equal to passed myPoint
     * @param b myPoint to check
     * @return whether they are equal
     */
    @Override
    public boolean equals(Object b){
        if (this == b) return true;
        if (b instanceof myPoint v) {return ((this.x == v.x) && (this.y == v.y) && (this.z == v.z));}
        return false;            
    }//equals
    
    /**
     * Linearly interpolate between two floating point values
     * @param a
     * @param s interpolant
     * @param b
     * @return
     */
    protected final static double _linInterp(double a, double s, double b) {return (1-s)*a + (s)*b;}
    
    /**
     * Linearly interpolate between two floating point values, restricting result to lie in [min, max]
     * @param a
     * @param s interpolant
     * @param b
     * @param min
     * @param max
     * @return
     */
    protected final static double _cappedLinInterp(double a, double s, double b, double min, double max) {
        double res = (1-s)*a + (s)*b; 
        return (res > max ? max : res < min ? min : res);
    }
    
    public String toStrCSV(){return toStrCSV("%.4f");}    
    public String toStrCSV(String fmt){return "" + String.format(fmt,this.x) + ", " + String.format(fmt,this.y) + ", " + String.format(fmt,this.z);}    
    public String toStrBrf(){return "(" + String.format("%.4f",this.x) + ", " + String.format("%.4f",this.y) + ", " + String.format("%.4f",this.z)+")";}    
    public String toString(){return "|(" + String.format("%.4f",this.x) + ", " + String.format("%.4f",this.y) + ", " + String.format("%.4f",this.z)+")";}
}
