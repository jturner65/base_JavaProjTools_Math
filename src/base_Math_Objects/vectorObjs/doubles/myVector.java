package base_Math_Objects.vectorObjs.doubles;

import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.floats.myVectorf;

public class myVector extends myPoint{
    public double sqMagn, magn;
    /**
     * vector constants available to all consumers of myVector
     */   
    /**
     * zero vector
     */
    public static final myVector ZEROVEC = new myVector(0,0,0);
    /**
     * Up for Graphics framing (with z as depth, and y increasing downward)
     */
    public static final myVector GRAPHICS_UP = new myVector(0,-1,0);
    /**
     * Right for Graphics framing (with z as depth, and y increasing downward)
     */
    public static final myVector GRAPHICS_RIGHT = new myVector(1,0,0);
    /**
     * Forward for Graphics framing (with z as depth, and y increasing downward
     */
    public static final myVector GRAPHICS_FORWARD = new myVector(0,0,1);
    /**
     * Physics frame up as z
     */
    public static final myVector UP = new myVector(0,0,1);
    /**
     * Physics frame right as y
     */
    public static final myVector RIGHT = new myVector(0,1,0);
    /**
     * Physics frame forward as x
     */
    public static final myVector FORWARD = new myVector(1,0,0);
    /**
     * 
     * @param _x
     * @param _y
     * @param _z
     */    
    public myVector(double _x, double _y, double _z){super(_x,_y,_z); this._mag();}  
    /**
     * Build a vector using the first 3 values in the passed array. Will fail if less than 3 values.
     * @param vals must be 3+ values in length. Any values past the first 3 will be ignored.
     */
    public myVector(float[] vals){super(vals);this._mag();}    
    /**
     * Build a vector using the first 3 values in the passed array. Will fail if less than 3 values.
     * @param vals must be 3+ values in length. Any values past the first 3 will be ignored.
     */
    public myVector(double [] vals){super(vals);this._mag();} 
    /**
     * Copy Constructor
     * @param p
     */
    public myVector(myVector p){super(p.x, p.y, p.z); this.magn=p.magn; this.sqMagn=p.sqMagn; }  
    /**
     * Copy Constructor
     * @param p
     */
    public myVector(myVectorf p){super(p.x, p.y, p.z); this.magn=p.magn; this.sqMagn=p.sqMagn; }  
    /**
     * Empty constructor
     */
    public myVector(){ }
    /**
     * Create a vector from point a to point b
     * @param a
     * @param b
     */
    public myVector(myPointf a, myPointf b){this(b.x-a.x,b.y-a.y,b.z-a.z);}
    /**
     * Create a vector from point a to point b
     * @param a
     * @param b
     */
    public myVector(myPoint a, myPointf b){this(b.x-a.x,b.y-a.y,b.z-a.z);}
    /**
     * Create a vector from point a to point b
     * @param a
     * @param b
     */
    public myVector(myPoint a, myPoint b){this(b.x-a.x,b.y-a.y,b.z-a.z);}
    
    /**
     * Create a vector from origin to point a
     * @param a
     */
    public myVector(myPoint a){this(a.x,a.y,a.z);}
    /**
     * Create a vector that interpolates between vector a and vector b by _s value
     * @param a
     * @param _s
     * @param b
     */
    public myVector(myVector a, double _s, myVector b) {super(a,_s,b);this._mag();    }
    /**
     * build unit vector copy of v
     * @param v
     * @return
     */
    public static myVector _unit(myVector v){myVector u = new myVector(v); return u._normalize(); }
    /**
     * Create unit vector that interpolates between vector a and vector b by _s value
     * @param a
     * @param _s
     * @param b
     * @return
     */
    public static myVector _unit(myVector a, double _s, myVector b){myVector r = new myVector(a,_s,b); return r._normalize(); }
    /**
     * build unit vector pointing from a to point b
     * @param a 
     * @param b
     * @return
     */
    public static myVector _unit(myPoint a, myPoint b){myVector u = new myVector(a,b); return u._normalize(); }
    /**
     * build unit vector from origin to v
     * @param v
     * @return
     */
    public static myVector _unitFromPoint(myPoint v){myVector u = new myVector(v); return u._normalize(); }
    /**
     * Clear the values in this vector
     */
    @Override
    public void clear() {super.clear();this.magn = 0; this.sqMagn=0;}
    /**
     * Set the values in this vector
     */
    @Override
    public void set(double _x, double _y, double _z){ super.set(_x, _y, _z); this._mag(); }
    /**
     * Set the values in this vector to be the passed vector's values
     */
    public void set(myVector p){ this.x = p.x; this.y = p.y; this.z = p.z;  this.magn=p.magn; this.sqMagn=p.sqMagn; }
    /**
     * Set the values in this vector to be the vector from point p to point q
     */
    public void set(myPoint p, myPoint q){ this.x = q.x - p.x; this.y = q.y - p.y; this.z = q.z - p.z;  this._mag();}
    /**
     * Set all the values in this vector - DOES NOT RECALCULATE MAGNITUDE
     * @param _x
     * @param _y
     * @param _z
     * @param _sqMagn
     */
    public void set(double _x, double _y, double _z, double _sqMagn){ super.set(_x, _y, _z); this.sqMagn = _sqMagn; }
    /**
     * Return average vector of this vector and passed vector
     * @param q
     * @return
     */
    public myVector _avgWithMe(myVector q) {return new myVector((this.x+q.x)/2.0,(this.y+q.y)/2.0,(this.z+q.z)/2.0);} 
    /**
     * Static Method : Return average of the two passed vectors
     * @param p
     * @param q
     * @return
     */
    public static myVector _average(myVector p, myVector q) {return new myVector((p.x+q.x)/2.0,(p.y+q.y)/2.0,(p.z+q.z)/2.0);} 
    /**
     * Multiply this vector by n. Returns this vector
     */
    @Override
    public myVector _mult(double n){ super._mult(n); this._mag(); return this; }  
    /**
     * Static Method : Multiply passed vector p by n and returns result
     * @param p
     * @param n
     * @return
     */
    public static myVector _mult(myVector p, double n){return new myVector(p.x * n, p.y * n, p.z * n); }
    /**
     * Static Method : Multiply passed vector p by n and returns result
     * @param p
     * @param n
     * @return
     */
    public static myVector _mult(myVectorf p, double n){ return new myVector(p.x * n, p.y * n, p.z * n); }
    /**
     * Static Method : Elementwise-multiply passed vector p by passed vector q and returns resultant vector
     * @param p
     * @param q
     * @return
     */
    public static myVector _mult(myVector p, myVector q){ return new myVector(p.x *q.x, p.y * q.y, p.z * q.z); }
    /**
     * Static Method : Elementwise-multiply passed vector p by passed vector q and returns resultant vector
     * @param p
     * @param q
     * @return
     */
    public static myVector _ewise_mult(myVector p, myVector q){ return new myVector(p.x *q.x, p.y * q.y, p.z * q.z); }
    /**
     * Static Method : Elementwise-multiply first two arguments and put result in 3rd argument
     * @param p
     * @param q
     * @param r
     */
    public static void _mult(myVector p, myVector q, myVector r){ myVector result = new myVector(p.x *q.x, p.y * q.y, p.z * q.z); r.set(result);}  
    /**
     * Divide this vector by q
     */
    @Override
    public void _div(double q){super._div(q); this._mag();}  
    /**
     * Static Method : Divide vector p by n
     * @param p
     * @param n
     * @return
     */
    public static myVector _div(myVector p, double n){ if(n==0) return p; return new myVector(p.x / n, p.y / n, p.z / n); }
    /**
     * Static Method : Element-wise divide vector p by vector q
     * @param p
     * @param q
     * @return
     */
    public static myVector _div(myVector p, myVector q){ if((q.x==0)||(q.y==0)|| (q.z==0))return p; return new myVector(p.x / q.x, p.y / q.y, p.z / q.z); }
    /**
     * Static Method : Element-wise divide vector p by vector q
     * @param p
     * @param q
     * @return
     */
    public static myVector _ewise_div(myVector p, myVector q){ if((q.x==0)||(q.y==0)|| (q.z==0))return p; return new myVector(p.x / q.x, p.y / q.y, p.z / q.z); }
    /**
     * Add coordinate values to this vector
     */
    @Override
    public void _add(double _x, double _y, double _z){ super._add(_x, _y, _z); this._mag(); }
    /**
     * Add a vector to this vector
     * @param v
     */
    public void _add(myVector v){ this.x += v.x; this.y += v.y; this.z += v.z;  this._mag();  }
    /**
     * Add a vector to this vector
     * @param v
     */
    public void _add(myVectorf v){ this.x += v.x; this.y += v.y; this.z += v.z;  this._mag();  }
    /**
     * Static Method : add two vectors and return result
     * @param p,q : vectors to add
     * @return resulting point of element-wise addition of p + q
     */
    public static myVector _add(myVector p, myVector q){ return new myVector(p.x + q.x, p.y + q.y, p.z + q.z); }                   
    /**
     * Static Method : add vector to passed point and return result
     * @param p point
     * @param q vector
     * @return resulting point of element-wise addition of p + q
     */
    public static myVector _add(myPoint p, myVector q){ return new myVector(p.x + q.x, p.y + q.y, p.z + q.z); }
    /**
     * Static Method : add two vectors, putting result in 3rd argument
     * @param p,q : vectors to add
     * @param r : resulting point of element-wise addition of p + q
     */
    public static void _add(myVector p, myVector q, myVector r){ myVector result = new myVector(p.x + q.x, p.y + q.y, p.z + q.z); r.set(result);}  
    /**
     * Static Method : Find the sum of an array of vectors, minimizes extra calcs.
     * @param vAra
     * @return
     */
    public static myVector _add(myVectorf[] vAra) {
 // aggregate to minimize magnitude calcs
        float x = 0, y = 0, z = 0;
        for(int i=0;i<vAra.length;++i) {
            x += vAra[i].x;
            y += vAra[i].y;
            z += vAra[i].z;
        }
        return new myVector(x,y,z);
    }    
    /**
     * Subtract values from this vector
     */
    @Override
    public void _sub(double _x, double _y, double _z){ super._sub(_x, _y, _z);  this._mag(); }
    /**
     * Element-wise subtract passed vector from this vector
     * @param v
     */
    public void _sub(myVector v){ this.x -= v.x; this.y -= v.y; this.z -= v.z;  this._mag(); }
    /**
     * Static Method : Subtract q from p and return result
     * @param p
     * @param q
     * @return
     */
    public static myVector _sub(myPoint p, myPoint q){ return new myVector(p.x - q.x, p.y - q.y, p.z - q.z);}
    /**
     * Static Method : Subtract q from p, putting result in 3rd argument
     * @param p 
     * @param q
     * @param r
     */
    public static void _sub(myVector p, myVector q, myVector r){ myVector result = new myVector(p.x - q.x, p.y - q.y, p.z - q.z); r.set(result);}  
    /**
     * Calc, set and return this vector's magnitude
     * @return
     */
    public double _mag(){ this.magn = Math.sqrt(this._SqMag()); return this.magn; }  
    /**
     * Calc, set and return this vector's sq magnitude
     * @return
     */
    public double _SqMag(){ this.sqMagn = ((this.x*this.x) + (this.y*this.y) + (this.z*this.z)); return this.sqMagn; }
    /**
     * Sets length of vector
     * @param _newMag
     */
    public void _scale(double _newMag){this._normalize()._mult(_newMag);}
    /**
     * Normalize this vector and return it
     * @return
     */
    public myVector _normalize(){this._mag();if(magn==0){return this;}this.x /= magn; this.y /=magn; this.z /=magn; _mag();return this;}
    /**
     * Static Method : build a normalizeed version of the passed vector
     * @param v
     * @return
     */
    public static myVector _normalize(myVector v){double magn = v._mag(); if(magn==0){return v;} return new myVector( v.x / magn, v.y / magn, v.z / magn); }
    /**
     * Recalculate and return this vector's magnitude
     * @return
     */
    public double _norm(){return _mag();}
    /**
     * Static Method : Calculate and return the L2 norm of the passed vector
     * @param v
     * @return
     */    
    public static double _L2Norm(myVector v){return Math.sqrt(v._SqMag());}
    /**
     * Static Method : Calculate and return the squared L2 norm of the passed vector
     * @param v
     * @return
     */
    public static double _L2SqNorm(myVector v){return v._SqMag();}
    /**
     * Build new normalized version of this vector (this vector remains unchanged)
     * @return
     */
    public myVector _normalized(){double magn = this._mag(); myVector newVec = (magn == 0) ? (new myVector(0,0,0)) : (new myVector( this.x / magn, this.y / magn, this.z / magn)); newVec._mag(); return newVec;}
    /**
     * Build a copy of this vector
     */
    @Override
    public myVector cloneMe(){myVector retVal = new myVector(this.x, this.y, this.z); return retVal;}  
    /**
     * Static Method : find the squared distance between the two passed vectors 
     * @param q
     * @param r
     * @return
     */ 
    public static double _SqrDist(myVector q, myVector r){  return ((r.x - q.x)*(r.x - q.x)) + ((r.y - q.y)*(r.y - q.y)) + ((r.z - q.z)*(r.z - q.z));}
    
    /**
     * Static Method : find the distance between the two passed vectors 
     * @param q
     * @param r
     * @return
     */
    public static double _dist(myVector q, myVector r){  return Math.sqrt(((r.x - q.x) *(r.x - q.x)) + ((r.y - q.y) *(r.y - q.y)) + ((r.z - q.z) *(r.z - q.z)));}
    /**
     * Static Method : find the distance between the passed vector and the vector represented by the passed 3 values
     * @param r
     * @param qx
     * @param qy
     * @param qz
     * @return
     */
    public static double _dist(myVector r, double qx, double qy, double qz){  return Math.sqrt(((r.x - qx) *(r.x - qx)) + ((r.y - qy) *(r.y - qy)) + ((r.z - qz) *(r.z - qz)));}    
    /**
     * Find the cross product of this vectoor and the passed vector
     * @param b
     * @return
     */        
    public myVector _cross(myVector b){ return new myVector((this.y * b.z) - (this.z*b.y), (this.z * b.x) - (this.x*b.z), (this.x * b.y) - (this.y*b.x));}
    /**
     * Static Method : Find the cross product of the two passed vectors
     * @param a
     * @param b
     * @return
     */
    public static myVector _cross(myVector a, myVector b){        return a._cross(b);}
    /**
     * Static Method : Find the cross product of the vectors formed by each set of 3 values
     * @param ax
     * @param ay
     * @param az
     * @param bx
     * @param by
     * @param bz
     * @return
     */   
    public static myVector _cross(double ax, double ay, double az, double bx, double by, double bz){        return new myVector((ay*bz)-(az*by), (az*bx)-(ax*bz), (ax*by)-(ay*bx));}
    /**
     * Find the dot product of this vector dotted with the passed vector
     * @param b
     * @return
     */   
    public double _dot(myVector b){return ((this.x * b.x) + (this.y * b.y) + (this.z * b.z));}
    /**
     * Find the dot product of this vector dotted with the passed vector
     * @param b
     * @return
     */
    public double _dot(myVectorf b){return ((this.x * b.x) + (this.y * b.y) + (this.z * b.z));}
    /**
     * Static Method : Find the dot product of two two passed vectors
     * @param a
     * @param b
     * @return
     */
    public static double _dot(myVector a, myVector b){        return a._dot(b);}
    /**
     * Static Method : U|V det product
     * @param U
     * @param V
     * @return
     */
    public static double _det3(myVector U, myVector V) {double udv = U._dot(V); return (Math.sqrt((U.sqMagn*V.sqMagn) - (udv*udv))); }
    /**
     * Static Method : U*(VxW)  mixed product, determinant - measures 6x the volume of the parallelapiped formed by passed vectors
     * @param U
     * @param V
     * @param W
     * @return
     */
    public static double _mixProd(myVector U, myVector V, myVector W) {return U._dot(myVector._cross(V,W)); }  
    /**
     * Static Method : U * (VxW)>0  U,V,W are clockwise
     * @param U
     * @param V
     * @param W
     * @return
     */
    public static boolean _isCW_Vecs(myVector U, myVector V, myVector W) {return _mixProd(U,V,W)>0; } 
    /**
     * Static Method : Area of triangle described by 3 points 
     * @param A, B, C Triangle verts
     * @return area of proscribed triangle
     */
    public static double _area(myPoint A, myPoint B, myPoint C) {    myVector x = new myVector(A,B), y = new myVector(A,C), z = x._cross(y);     return z.magn/2.0; } 
    /**
     * Static Method : Returns volume of tetrahedron defined by A,B,C,D
     * @param A,B,C,D verts of tet
     * @return volume
     */
    public static double _volume(myPoint A, myPoint B, myPoint C, myPoint D) {return _mixProd(new myVector(A,B),new myVector(A,C),new myVector(A,D))/6.0; }
    /**
     * Static Method : Returns true if tet is oriented so that A sees B, C, D clockwise
     * @param A,B,C,D verts of tet
     * @return if tet is oriented clockwise (A->B->C->D)
     */
    public static boolean _isCW_Tet(myPoint A, myPoint B, myPoint C, myPoint D) {return _volume(A,B,C,D)>0; } 
    /**
     * Static Method : Returns true if U and V are almost parallel
     * @param U
     * @param V
     * @return
     */
    public static boolean isParallel(myVector U, myVector V) {return U._cross(V).magn<U.magn*V.magn*MyMathUtils.EPS; }
    
    /**
     * Static Method : This will return a unit normal ortho to the plane spanned by U and V. If U and V are parallel, will provide a normal to U.
     * @param U 
     * @param V
     * @return
     */
    public static myVector _findNormToPlane(myVector U, myVector V) {
        myVector norm = V._cross(U);
        if(norm.magn<U.magn*V.magn*MyMathUtils.EPS) {
 //parallel, find normal to U and modified V 
 //that will never be coincident no matter what V is
            norm = myVector._cross(U.x, U.y, U.z, (V.x*2)+1, (V.y*3)+2, (V.z*4)+3);            
        }
        norm._normalize();
        return norm;
    }
    /**
     * Static Method : Returns the angle between the two passed vectors
     * @param v1
     * @param v2
     * @return
     */
    public static double _angleBetween(myVector v1, myVector v2) {
        double _v1Mag = v1._mag(), 
                _v2Mag = v2._mag(), 
                dotProd = v1._dot(v2),
                cosAngle = dotProd/(_v1Mag * _v2Mag),
                angle = Math.acos(cosAngle);
        return angle;
    }//_angleBetween
    /**
     * Returns the angle between this and the passed vector
     * @param v2
     * @return
     */   
    public double angleWithMe(myVector v2) {
        double _v1Mag = _mag(), 
                _v2Mag = v2._mag(), 
                dotProd = _dot(v2),
                cosAngle = dotProd/(_v1Mag * _v2Mag),
                angle = Math.acos(cosAngle);
        return angle;
    }//_angleBetween    
    /**
     * Static Method : Rotate v1 around axis unit vector u, by Pi/2, around origin
     * @param v1
     * @param u
     * @return
     */
    public static myVector _rotAroundAxis(myVector v1, myVector u){return _rotAroundAxis(v1, u, MyMathUtils.HALF_PI);}
    /**
     * Static Method : Rotate v1 around axis unit vector u, by given angle thet, around origin
     * @param v1
     * @param u
     * @param thet
     * @return
     */
    public static myVector _rotAroundAxis(myVector v1, myVector u, double thet){        
        double cThet = Math.cos(thet), sThet = Math.sin(thet), oneMC = 1-cThet,
                ux2 = u.x*u.x, uy2 = u.y*u.y, uz2 = u.z*u.z,
                uxy = u.x * u.y, uxz = u.x * u.z, uyz = u.y*u.z,
                uzS = u.z*sThet, uyS = u.y*sThet, uxS = u.x*sThet,
                uxzC1 = uxz *oneMC, uxyC1 = uxy*oneMC, uyzC1 = uyz*oneMC;
        //build rot matrix in vector def
        myVector res = new myVector(
                (ux2*oneMC+cThet) * v1.x + (uxyC1-uzS)         * v1.y + (uxzC1+uyS) *v1.z,
                (uxyC1+uzS)       * v1.x + (uy2*oneMC+cThet)* v1.y + (uyzC1-uxS) *v1.z,
                (uxzC1-uyS)       * v1.x + (uyzC1+uxS)        * v1.y + (uz2*oneMC + cThet) * v1.z);
        
        return res;        
    }
    /**
     * Rotate this vector around given axis vector u by angle thet
     * @param u
     * @param thet
     * @return
     */
    public myVector rotMeAroundAxis(myVector u, double thet){        
        double cThet = Math.cos(thet), sThet = Math.sin(thet), oneMC = 1-cThet,
                ux2 = u.x*u.x, uy2 = u.y*u.y, uz2 = u.z*u.z,
                uxy = u.x * u.y, uxz = u.x * u.z, uyz = u.y*u.z,
                uzS = u.z*sThet, uyS = u.y*sThet, uxS = u.x*sThet,
                uxzC1 = uxz *oneMC, uxyC1 = uxy*oneMC, uyzC1 = uyz*oneMC;
        //build rot matrix in vector def
        myVector res = new myVector(
                (ux2*oneMC+cThet) * this.x + (uxyC1-uzS)         * this.y + (uxzC1+uyS) *this.z,
                (uxyC1+uzS)       * this.x + (uy2*oneMC+cThet)* this.y + (uyzC1-uxS) *this.z,
                (uxzC1-uyS)       * this.x + (uyzC1+uxS)        * this.y + (uz2*oneMC + cThet) * this.z);
        
        return res;        
    }    

    /**
     * undefined for a vector - use rotMeAroundAxis
     */
    public final myPoint rotMeAroundPt(float a, myVector I, myVector J, myPoint G) {return this;}
    /**
     * undefined for a vector - rotMeAroundAxis
     */
    public final myPoint rotMeAroundPt(myPoint C, myPoint P, myPoint R)  {return this;}

    /**
     * Static Method : Find the angle between two vectors - Note this version is for 2D - relies on neither vector being coplanar with (0,0,1);
     * @param U
     * @param V
     * @return
     */    
    public static double _angleBetween_Xprod(myVector U, myVector V){
        return _angleBetween_Xprod(U,V, new myVector(0,0,1));
    }

    /**
     * Static Method : Find the angle between two vectors, with axis about which to determine sign
     * @param U
     * @param V
     * @param axis unit length axis about which to determine sign
     * @return
     */    
    public static double _angleBetween_Xprod(myVector U, myVector V, myVector axis){
        myVector cross = U._cross(V);
        double dot = U._dot(V);
        
        double angle = (float) Math.atan2(cross.magn,dot),
                sign = _mixProd(U,V,axis);
        if(sign<0){    angle=-angle;}    
        return angle;
    }

    /**
     * Returns if this vector is equal to passed vector
     * @param b vector to check
     * @return whether they are equal
     */
    @Override
    public boolean equals(Object b){
        if (this == b) return true;
        if (b instanceof myVector v) {return ((this.x == v.x) && (this.y == v.y) && (this.z == v.z));}
        return false;
    }
    /**
     * 
     */
    public String toStrCSV(){return toStrCSV("%.4f");}
    /**
     * 
     */
    public String toStrCSV(String fmt){return super.toStrCSV(fmt) + ", " + String.format(fmt,this.magn) + ", " + String.format(fmt,this.sqMagn);}
    /**
     * 
     */
    public String toStrBrf(){return super.toStrBrf() + ", " + String.format("%.4f",this.magn) + ", " + String.format("%.4f",this.sqMagn);}    
    /**
     * 
     */
    public String toString(){return super.toString()+ " | Mag:" + String.format("%.4f",this.magn)+ " | sqMag:" + String.format("%.4f",this.sqMagn);}
}// myVector