package base_Math_Objects.vectorObjs.floats;

import base_Math_Objects.vectorObjs.doubles.myPoint;

public class myPointf {
    /**
     * this point object's x coord
     */
    public float x;
    /**
     * this point object's y coord
     */
    public float y;
    /**
     * this point object's z coord
     */
    public float z;
    /**
     * Static final member : origin point in 3D
     */
    public static final myPointf ZEROPT = new myPointf(0,0,0);

    /**
     * build a point with given coordinates
     * @param _x : x coord
     * @param _y : y coord
     * @param _z : z coord
     */
    public myPointf(float _x, float _y, float _z){this.x = _x; this.y = _y; this.z = _z;}         //constructor 3 args  
    /**
     * build a point with given coordinates(as doubles)
     * @param _x : x coord
     * @param _y : y coord
     * @param _z : z coord
     */
    public myPointf(double _x, double _y, double _z){this.x = (float) _x; this.y = (float) _y; this.z = (float) _z;}         //constructor 3 args  
    /**
     * Build a point using the first 3 values in the passed array. Will fail if less than 3 values.
     * @param vals must be 3+ values in length. Any values past the first 3 will be ignored.
     */
    public myPointf(float[] vals){this.x = vals[0]; this.y = vals[1]; this.z = vals[2];}    
    /**
     * Build a point using the first 3 values in the passed array. Will fail if less than 3 values.
     * @param vals must be 3+ values in length. Any values past the first 3 will be ignored.
     */
    public myPointf(double[] vals){this.x = (float)vals[0]; this.y = (float)vals[1]; this.z = (float)vals[2];} 
    /**
     * copy constructor
     * @param p : point object to copy
     */
    public myPointf(myPointf p){ this(p.x, p.y, p.z); }                                                                                                               //constructor 1 arg  
    /**
     * build point as displacement from point A by vector B
     * @param A : starting point
     * @param B : displacement vector
     */
    public myPointf(myPointf A, myVectorf B) {this(A.x+B.x,A.y+B.y,A.z+B.z); };
    /**
     * Interpolate between A and B by s -> (0->1)
     * @param A : first point to interpolate from
     * @param s : value [0,1] to determine linear interpolation
     * @param B : second point to interpolate from
     */
    public myPointf(myPointf A, float s, myPointf B) {  this(_linInterp(A.x, s, B.x), _linInterp(A.y, s, B.y), _linInterp(A.z, s, B.z)); };        //builds a point somewhere in between a and b
    /**
     * empty constructor
     */
    public myPointf(){ this(0,0,0);}                                                                                                                               //constructor 0 args
    /**
     * clear this point's set values
     */
    public void clear() {this.x = 0; this.y = 0; this.z = 0;}    
    /**
     * Set this object's coordinate values
     * @param _x, _y, _z : new x,y,z coords of this object
     */
    public void set(float _x, float _y, float _z){ this.x = _x;  this.y = _y;  this.z = _z; }                                               //set 3 args 
    /**
     * Set this object's coordinate values (converting from doubles)
     * @param _x, _y, _z : new x,y,z coords of this object
     */
    public void set(double _x, double _y, double _z){ this.set((float)_x,(float)_y,(float)_z); }                                               //set 3 args 
    /**
     * Set this point's values as a copy of the passed point
     * @param p : the point to copy
     * @return this point
     */
    public myPointf set(myPointf p){ this.x = p.x; this.y = p.y; this.z = p.z; return this;}                                                                   //set 1 args
    /**
     * build and return the average (midpoint) of this point and the passed point
     * @param q : the point to find the midpoint with
     * @return the midpoint between this and q
     */    
    public myPointf _avgWithMe(myPointf q) {return new myPointf((this.x+q.x)/2.0f,(this.y+q.y)/2.0f,(this.z+q.z)/2.0f);} 
    /**
     * build and return the average (midpoint) of this point and the passed points
     * @param q,r : the points to find the midpoint with
     * @return the midpoint between this and q and r
     */    
    public myPointf _avgWithMe(myPointf q, myPointf r) {return new myPointf((this.x+q.x+r.x)/3.0f,(this.y+q.y+r.y)/3.0f,(this.z+q.z+r.z)/3.0f);} 
    /**
     * build and return the average (midpoint) of this point and the passed points
     * @param q,r : the points to find the midpoint with
     * @return the midpoint between this and q and r
     */    
    public myPointf _avgWithMe(myPointf q, myPointf r, myPointf s) {return new myPointf((this.x+q.x+r.x+s.x)/4.0f,(this.y+q.y+r.y+s.y)/4.0f,(this.z+q.z+r.z+s.z)/4.0f);} 
    /**
     * Static method : build the midpoint between the two passed points
     * @param p,q : 2 points to find the midpoint between
     * @return midpoint between p and q
     */
    public static myPointf _average(myPointf p, myPointf q) {return new myPointf((p.x+q.x)/2.0f,(p.y+q.y)/2.0f,(p.z+q.z)/2.0f);} 
    /**
     * Static method : build the midpoint between the three passed points
     * @param p,q,r : 3 points to find the midpoint between
     * @return midpoint between  p,q,r
     */
    public static myPointf _average(myPointf p, myPointf q, myPointf r) {return new myPointf((p.x+q.x+r.x)/3.0f,(p.y+q.y+r.y)/3.0f,(p.z+q.z+r.z)/3.0f);} 
    /**
     * Static method : build the midpoint between the two passed points
     * @param p,q,r,s : 4 points to find the midpoint between
     * @return midpoint between p,q,r,s
     */
    public static myPointf _average(myPointf p, myPointf q, myPointf r, myPointf s) {return new myPointf((p.x+q.x+r.x+s.x)/4.0f,(p.y+q.y+r.y+s.y)/4.0f,(p.z+q.z+r.z+s.z)/4.0f);} 
    /**
     * Static method : build the center of verticies for the array of passed points (COV)
     * @param pts : array of points
     * @return COV for pts array
     */
    public static myPointf _average(myPointf[] pts) {myPointf res = new myPointf();    for(myPointf p : pts) {    res._add(p);}return myPointf._div(res, 1.0f*pts.length);} 
    /**
     * multiply this point by n, modifying this point to be result
     * @param n : value to scale this vector by
     * @return this point, after scaling
     */
    public myPointf _mult(float n){ this.x *= n; this.y *= n; this.z *= n; return this; }                                                     //_mult 3 args  
    /**
     * Static method : return result of multiplying a point by a scalar
     * @param p : point to be scaled
     * @param n : scale value (to multiply)
     * @return point result of element-wise multiplication of p by n
     */
    public static myPointf _mult(myPointf p, float n){ return new myPointf(p.x * n, p.y * n, p.z * n);}                          //1 pt, 1 double
    /**
     * Static method : element-wise multiplication of two points, returning result
     * @param p,q : two points to multiply element-wise with each other
     * @return : element-wise product of p and q : (p.x * q.x) i + (p.y*q.y) j + (p.z*q.z) k
     */
    public static myPointf _mult(myPointf p, myPointf q){ return new myPointf(p.x *q.x, p.y * q.y, p.z * q.z);}           //return elementwise product
    /**
     * Static method : element-wise multiplication of two points, returning result
     * @param p,q : two points to element-wise multiply
     * @param r : destination for result of element-wise multiplication
     */
    public static void _mult(myPointf p, myPointf q, myPointf r){ myPointf result = new myPointf(p.x *q.x, p.y * q.y, p.z * q.z); r.set(result);}           //2 pt src, 1 pt dest      
    /**
     * divide this point by q, making this point equal to result.  No value checking is performed
     * @param q scalar value to divide this point by
     */
    public void _div(float q){this.x /= q; this.y /= q; this.z /= q; }  
    /**
     * Static Method : divide passed point p by n, making this point equal to result.  No value checking is performed
     * @param p : point to scale (divide)
     * @param n : value to divide p by
     * @return point result of per-element division of p by n
     */
    public static myPointf _div(myPointf p, float n){ if(n==0) return p; return new myPointf(p.x / n, p.y / n, p.z / n);}                          //1 pt, 1 double
    /**
     * add passed values to this point, making this point equal to result. 
     * @param _x : x coord to add to this.x
     * @param _y : y coord to add to this.y
     * @param _z : z coord to add to this.z
     */
    public void _add(float _x, float _y, float _z){ this.x += _x; this.y += _y; this.z += _z;   }                                            //_add 3 args
    /**
     * element-wise add passed point to this point, making this point equal to result. 
     * @param v point to add to this
     */
    public void _add(myPointf v){ this.x += v.x; this.y += v.y; this.z += v.z;   }  
    /**
     * Static Method : Add vector I to point O, returning result
     * @param O : origin point
     * @param I : displacement vector
     * @return : O + I
     */
    public static myPointf _add(myPointf O, myVectorf I){                                                        return new myPointf(O.x+I.x,O.y+I.y,O.z+I.z);}  
    /**
     * Static Method : Add vector I (scaled by a) to point O, returning result
     * @param O : origin point
     * @param a : scaling of I
     * @param I : displacement vector
     * @return : O + aI
     */
    public static myPointf _add(myPointf O, float a, myVectorf I){                                                return new myPointf(O.x+a*I.x,O.y+a*I.y,O.z+a*I.z);}                                        //2 vec
    /**
     * Static Method : Add vector I (scaled by a) and vector J (scaled by b) to point O, returning result
     * @param O : origin point
     * @param a : scaling of I
     * @param I : displacement vector
     * @param b : scaling of J
     * @param J : displacement vector
     * @return : O + aI + bJ
     */
    public static myPointf _add(myPointf O, float a, myVectorf I, float b, myVectorf J) {                        return new myPointf(O.x+a*I.x+b*J.x,O.y+a*I.y+b*J.y,O.z+a*I.z+b*J.z);}                      // O+xI+yJ
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
    public static myPointf _add(myPointf O, float a, myVectorf I, float b, myVectorf J, float c, myVectorf K) {    return new myPointf(O.x+a*I.x+b*J.x+c*K.x,O.y+a*I.y+b*J.y+c*K.y,O.z+a*I.z+b*J.z+c*K.z);} // O+xI+yJ+kZ
    /**
     * Static Method : add two points and return result
     * @param p,q : points to add
     * @return resulting point of element-wise addition of p + q
     */
    public static myPointf _add(myPointf p, myPointf q){return new myPointf(p.x + q.x, p.y + q.y, p.z + q.z); }
    /**
     * Static Method : add two points (one being a double-based point), putting result in 3rd argument
     * @param p : double-based point to add
     * @param q : float-based point to add
     * @return resulting point of element-wise addition of p + q
     */
    public static myPointf _add(myPoint p, myPointf q){return new myPointf(p.x + q.x, p.y + q.y, p.z + q.z); }
    /**
     * Static Method : add two points and return result
     * @param p,q : points to add
     * @param r : resulting point of element-wise addition of p + q
     */
    public static void _add(myPointf p, myPointf q, myPointf r){ myPointf result = new myPointf(p.x + q.x, p.y + q.y, p.z + q.z); r.set(result);}           //2 pt src, 1 pt dest  
    /**
     * subtract passed values to this point, making this point equal to result. 
     * @param _x : x coord to subtract from this.x
     * @param _y : y coord to subtract from this.y
     * @param _z : z coord to subtract from this.z
     */
    public void _sub(float _x, float _y, float _z){ this.x -= _x; this.y -= _y; this.z -= _z;  }                                                                   //_sub 3 args
    /**
     * element-wise subtract passed point from this point, making this point equal to result. 
     * @param v point to subtract this
     */    
    public void _sub(myPointf v){ this.x -= v.x; this.y -= v.y; this.z -= v.z;  }                                                                           //_sub 1 arg 
    /**
     * Static Method : subtract two points and return result
     * @param p,q : points to subtract
     * @return resulting point of element-wise subtraction of p - q
     */
    public static myPointf _sub(myPointf p, myPointf q){ return new myPointf(p.x - q.x, p.y - q.y, p.z - q.z);}
    /**
     * Static Method : subtract two points, putting result in 3rd argument
     * @param p,q : points to subtract
     * @param r : resulting point of element-wise subtraction of p - q
     */
    public static void _sub(myPointf p, myPointf q, myPointf r){ myPointf result = new myPointf(p.x - q.x, p.y - q.y, p.z - q.z); r.set(result);}       //2 pt src, 1 pt dest      
    /**
     * create a new copy point of this point
     * @return
     */
    public myPointf cloneMe(){return new myPointf(this.x, this.y, this.z); }      

    /**
     * calculate L1 (Manhattan) distance between this point and passed point.  Manhattan distance is
     * Math.abs((this.x - q.x)) + Math.abs((this.y - q.y)) + Math.abs((this.z - q.z))
     * @param q : point to find distance from
     * @return : L1 distance from this point to q 
     */    
    public float _L1Dist(myPointf q){return Math.abs((this.x - q.x)) + Math.abs((this.y - q.y)) + Math.abs((this.z - q.z)); }
    /**
     * Static Method : calculate L1 (Manhattan) distance between passed points.  Manhattan distance is
     * Math.abs((r.x - q.x)) + Math.abs((r.y - q.y)) + Math.abs((r.z - q.z))
     * @param q,r : points to find distance between
     * @return : L1 distance from q to r
     */    
    public static float _L1Dist(myPointf q, myPointf r){ return q._L1Dist(r);}
    
    /**
     * find squared L2 (Euclidean) distance from this point to q.  Squared L2 distance is
     * ((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z)) 
     * @param q point to find distance from
     * @return squared L2 Distance from this point to q
     */
    public float _SqrDist(myPointf q){ float dx=(this.x-q.x), dy=(this.y-q.y), dz=(this.z-q.z);return ((dx*dx) + (dy*dy) + (dz*dz)); }
    /**
     * Static Method : find squared L2 (Euclidean) distance from point q to point r.  Squared L2 distance is
     * ((r.x - q.x)*(r.x - q.x)) + ((r.y - q.y)*(r.y - q.y)) + ((r.z - q.z)*(r.z - q.z)) 
     * @param q,r : points to find distance between
     * @return squared L2 Distance from q to r
     */
    public static float _SqrDist(myPointf q, myPointf r){ float dx=(r.x-q.x), dy=(r.y-q.y), dz=(r.z-q.z);return ((dx*dx) + (dy*dy) + (dz*dz)); }
    /**
     * find L2 (Euclidean) distance from this point to q.  L2 (Euclidean) distance is
     * sqrt(((r.x - q.x)*(r.x - q.x)) + ((r.y - q.y)*(r.y - q.y)) + ((r.z - q.z)*(r.z - q.z))) 
     * @param q point to find distance to
     * @return L2 Distance from this point to q
     */
    public float _dist(myPointf q){ float dx=(this.x-q.x), dy=(this.y-q.y), dz=(this.z-q.z);return (float) Math.sqrt(((dx*dx) + (dy*dy) + (dz*dz)));}
    /**
     * Static Method : find L2 (Euclidean) distance from point q to point r.  Squared L2 distance is
     * sqrt(((r.x - q.x)*(r.x - q.x)) + ((r.y - q.y)*(r.y - q.y)) + ((r.z - q.z)*(r.z - q.z))) 
     * @param q,r : points to find distance between
     * @return L2 Distance from q to r
     */
    public static float _dist(myPointf q, myPointf r){double dx=(r.x-q.x), dy=(r.y-q.y), dz=(r.z-q.z);return (float) Math.sqrt(((dx*dx) + (dy*dy) + (dz*dz)));}
    /**
     * find L2 (Euclidean) distance from this point to passed coordinates.  Squared L2 distance is
     * sqrt(((this.x - qx)*(this.x - qx)) + ((this.y - qy)*(this.y - qy)) + ((this.z - qz)*(this.z - qz)))
     * @param qx,qy,qz : coordinates to find distance to
     * @return L2 Distance from this to [qx,qy,qz]
     */
    public float _dist(float qx, float qy, float qz){ float dx=(this.x-qx), dy=(this.y-qy), dz=(this.z-qz);return (float) Math.sqrt(((dx*dx) + (dy*dy) + (dz*dz)));}
    /**
     * Static Method : find L2 (Euclidean) distance from point q to passed coordinates.  Squared L2 distance is
     * sqrt(((r.x - qx)*(r.x - qx)) + ((r.y - qy)*(r.y - qy)) + ((r.z - qz)*(r.z - qz)))
     * @param r : point to find distance from
     * @param qx,qy,qz : coordinates to find distance to
     * @return L2 Distance from r to [qx,qy,qz]
     */
    public static float _dist(myPointf r, float qx, float qy, float qz){ float dx=(r.x-qx), dy=(r.y-qy), dz=(r.z-qz);return (float) Math.sqrt((dx*dx) + (dy*dy) + (dz*dz)); }      
    /**
     * return the values of this point as an array of floats
     * @return array of floats {x,y,z}
     */
    public float[] asArray(){return new float[]{x,y,z};}
    /**
     * return the values of this point as an array of doubles
     * @return array of doubles {x,y,z}
     */
    public double[] asDblArray(){return new double[]{x,y,z};}
    /**
     * return the values of this point as a homogeneous point array
     * @return array of floats {x,y,z, 1}
     */
    public float[] asHAraPt(){return new float[]{this.x, this.y, this.z,1};}
    /**
     * return the values of this point as a homogeneous vector array
     * @return array of doubles {x,y,z, 0}
     */
    public float[] asHAraVec(){return new float[]{this.x, this.y, this.z,0};}
    
    /**
     * return the values of this point as a homogeneous point array
     * @return array of floats {x,y,z, 1}
     */
    public double[] asHAraPt_Dbl(){return new double[]{this.x, this.y, this.z,1};}
    /**
     * return the values of this point as a homogenous vector array
     * @return array of doubles {x,y,z, 0}
     */
    public double[] asHAraVec_Dbl(){return new double[]{this.x, this.y, this.z,0};}
    
    /**
     * this will calculate normalized barycentric coordinates (u, v, w) for pt w/respect to first 3 control points
     * @param cntlPts control points of poly
     * @return
     */
    public final float[] calcNormBaryCoords(myPointf[] cntlPts) {        
        //pt = u * cntlPts[0] + v * cntlPts[1] + w * cntlPts[2]
        myVectorf AB = new myVectorf(cntlPts[0],cntlPts[1]),
                AC = new myVectorf(cntlPts[0],cntlPts[2]),
                AP = new myVectorf(cntlPts[0],this);
        float d00 = AB.sqMagn, d01 = AB._dot(AC), 
            d11 = AC.sqMagn, d20 = AP._dot(AB), d21 = AP._dot(AC);
        float 
            denom = d00 * d11 - d01 * d01,
            v = (d11 * d20 - d01 * d21) / denom,
            w =  (d00 * d21 - d01 * d20) / denom,
            u = 1.0f - v- w;        
        return new float[] {u,v,w};
    }
    
    /**
     * given control points and passed normalized barycentric coordinates, calculate resultant point
     * @param cntlPts
     * @param pointNBC
     * @return
     */    
    public static final myPointf calcPointFromNormBaryCoords(myPointf[] cntlPts, float[] pointNBC) {
        myPointf res = new myPointf();
        for(int i=0;i<pointNBC.length;++i) {    res._add(myPointf._mult(cntlPts[i], pointNBC[i]));}    
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
    public myPointf rotMeAroundPt(float a, myVectorf I, myVectorf J, myPointf G) {
        double x= myVectorf._dot(new myVectorf(G,this),myVectorf._unit(I)), y=myVectorf._dot(new myVectorf(G,this),myVectorf._unit(J)); 
        double c=Math.cos(a), s=Math.sin(a); 
        float iXVal = (float) (x*c-x-y*s), jYVal= (float) (x*s+y*c-y);            
        return myPointf._add(this,iXVal,I,jYVal,J); 
    }
    
    /**
     * returns rotated version of this point by angle(CP,CR) parallel to plane (C,P,R)
     * @param C
     * @param P
     * @param R
     * @return this point rotated
     */
    public myPointf rotMeAroundPt(myPointf C, myPointf P, myPointf R) { // returns rotated version of Q by angle(CP,CR) parallel to plane (C,P,R)
        myVectorf I0=myVectorf._unit(C,P), I1=myVectorf._unit(C,R), V=new myVectorf(C,this); 
        double c=myPointf._dist(I0,I1), s=Math.sqrt(1.-(c*c)); 
        if(Math.abs(s)<0.00001) return this;        
        myVectorf J0=myVectorf._add(myVectorf._mult(I1,1./s),myVectorf._mult(I0,-c/s));  
        myVectorf J1=myVectorf._add(myVectorf._mult(I0,-s),myVectorf._mult(J0,c));  
        float x=V._dot(I0), y=V._dot(J0);  
        return myPointf._add(this,x,myVectorf._sub(I1,I0),y,myVectorf._sub(J1,J0)); 
    }
    
    /**
     * Returns whether the passed point is within eps of this point
     * @param p
     * @param eps
     * @return
     */
    public boolean clickIn(myPointf p, float eps) { return(_dist(p) < eps);}
    
    /**
     * returns if this pttor is equal to passed pttor
     * @param b myPointf to check
     * @return whether they are equal
     */
    public boolean equals(Object b){
        if (this == b) return true;
        if (b instanceof myPointf v) {return ((this.x == v.x) && (this.y == v.y) && (this.z == v.z));}
        return false;
    }
    
    /**
     * Linearly interpolate between two floating point values
     * @param a
     * @param s interpolant
     * @param b
     * @return
     */
    protected final static float _linInterp(float a, float s, float b) {return (1-s)*a + (s)*b;}
    
    /**
     * Linearly interpolate between two floating point values, restricting result to lie in [min, max]
     * @param a
     * @param s interpolant
     * @param b
     * @param min
     * @param max
     * @return
     */
    protected final static float _cappedLinInterp(float a, float s, float b, float min, float max) {
        float res = (1-s)*a + (s)*b; 
        return (res > max ? max : res < min ? min : res);
    }
   
    public String toStrCSV(){return toStrCSV("%.4f");}    
    public String toStrCSV(String fmt){return "" + String.format(fmt,this.x) + ", " + String.format(fmt,this.y) + ", " + String.format(fmt,this.z);}    
    public String toStrBrf(){return "(" + String.format("%.4f",this.x) + ", " + String.format("%.4f",this.y) + ", " + String.format("%.4f",this.z)+")";}    
    public String toString(){return "|(" + String.format("%.4f",this.x) + ", " + String.format("%.4f",this.y) + ", " + String.format("%.4f",this.z)+")";}
}