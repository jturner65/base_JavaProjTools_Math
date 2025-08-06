package base_Math_Objects.vectorObjs.floats;

import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.doubles.myPoint;
import base_Math_Objects.vectorObjs.doubles.myVector;

public class myVectorf extends myPointf{
    public float sqMagn, magn;
    
    /**
     * vector constants available to all consumers of myVectorf
     */
    /**
     * zero vector
     */
    public static final myVectorf ZEROVEC = new myVectorf(0,0,0);
    /**
     * Up for Graphics framing (with z as depth, and y increasing downward)
     */
    public static final myVectorf GRAPHICS_UP = new myVectorf(0,-1,0);
    /**
     * Right for Graphics framing (with z as depth, and y increasing downward)
     */
    public static final myVectorf GRAPHICS_RIGHT = new myVectorf(1,0,0);
    /**
     * Forward for Graphics framing (with z as depth, and y increasing downward
     */
    public static final myVectorf GRAPHICS_FORWARD = new myVectorf(0,0,1);
    /**
     * Physics frame up as z
     */
    public static final myVectorf UP = new myVectorf(0,0,1);
    /**
     * Physics frame right as y
     */
    public static final myVectorf RIGHT = new myVectorf(0,1,0);
    /**
     * Physics frame forward as x
     */
    public static final myVectorf FORWARD = new myVectorf(1,0,0);
    /**
     * 
     * @param _x
     * @param _y
     * @param _z
     */
    public myVectorf(float _x, float _y, float _z){super(_x,_y,_z); this._mag();}
    /**
     * 
     * @param _x
     * @param _y
     * @param _z
     */
    public myVectorf(double _x, double _y, double _z){super((float)_x,(float)_y,(float)_z); this._mag();}
    
    /**
     * Build a vector using the first 3 values in the passed array. Will fail if less than 3 values.
     * @param vals must be 3+ values in length. Any values past the first 3 will be ignored.
     */
    public myVectorf(float[] vals){super(vals);this._mag();}    
    /**
     * Build a vector using the first 3 values in the passed array. Will fail if less than 3 values.
     * @param vals must be 3+ values in length. Any values past the first 3 will be ignored.
     */
    public myVectorf(double[] vals){super(vals);this._mag();} 
    
    /**
     * Copy constructor
     * @param p
     */
    public myVectorf(myVectorf p){super(p.x, p.y, p.z); this.magn = p.magn;this.sqMagn = p.sqMagn;}
    /**
     * Copy Constructor
     * @param p
     */
    public myVectorf(myVector p){super((float)p.x, (float)p.y, (float)p.z); this.magn = (float) p.magn;this.sqMagn = (float) p.sqMagn;}  
    /**
     * Empty constructor
     */
    public myVectorf(){ }
    /**
     * Create a vector from point a to point b
     * @param a
     * @param b
     */
    public myVectorf(myPointf a, myPointf b){this(b.x-a.x,b.y-a.y,b.z-a.z);}
    /**
     * Create a vector from point a to point b
     * @param a
     * @param b
     */
    public myVectorf(myPoint a, myPointf b){this(b.x-a.x,b.y-a.y,b.z-a.z);}
    /**
     * Create a vector from point a to point b
     * @param a
     * @param b
     */
    public myVectorf(myPoint a, myPoint b){this(b.x-a.x,b.y-a.y,b.z-a.z);}
    
    /**
     * Create a vector from origin to point a
     * @param a
     */
    public myVectorf(myPointf a){this(a.x,a.y,a.z);}    
    /**
     * Create a vector that interpolates between vector a and vector b by _s value
     * @param a
     * @param _s
     * @param b
     */
    public myVectorf(myVectorf a, float _s, myVectorf b) {super(a,_s,b);this._mag();    }
    /**
     * build unit vector copy of v
     * @param v
     * @return
     */
    public static myVectorf _unit(myVectorf v){myVectorf u = new myVectorf(v); return u._normalize(); }
    /**
     * Create unit vector that interpolates between vector a and vector b by _s value
     * @param a
     * @param _s
     * @param b
     * @return
     */
    public static myVectorf _unit(myVectorf a, float _s, myVectorf b){myVectorf r = new myVectorf(a,_s,b); return r._normalize(); }
    /**
     * build unit vector pointing from a to point b
     * @param a 
     * @param b
     * @return
     */
    public static myVectorf _unit(myPointf a, myPointf b){myVectorf u = new myVectorf(a,b); return u._normalize(); }    
   /**
    * build unit vector from origin to v
    * @param v
    * @return
    */
    public static myVectorf _unitFromPoint(myPointf v){myVectorf u = new myVectorf(v); return u._normalize(); }
    /**
     * Clear the values in this vector
     */
    @Override
    public void clear() {super.clear();this.magn = 0; this.sqMagn=0;}
    /**
     * Set the values in this vector
     */
    @Override
    public void set(float _x, float _y, float _z){ super.set(_x, _y, _z); this._mag(); }
    /**
     * Set the values in this vector
     */
    @Override
    public void set(double _x, double _y, double _z){ this.set((float)_x,(float)_y,(float)_z); }//set 3 args
    /**
     * Set the values in this vector to be the passed vector's values
     */
    public void set(myVectorf p){ this.x = p.x; this.y = p.y; this.z = p.z; this.magn=p.magn; this.sqMagn=p.sqMagn; }
    /**
     * Set the values in this vector to be the vector from point p to point q
     */
    public void set(myPointf p, myPointf q){ this.x = q.x - p.x; this.y = q.y - p.y; this.z = q.z - p.z;  this._mag();}
    /**
     * Set all the values in this vector - DOES NOT RECALCULATE MAGNITUDE
     * @param _x
     * @param _y
     * @param _z
     * @param _sqMagn
     */
    public void set(float _x, float _y, float _z, float _sqMagn){ super.set(_x, _y, _z); this.sqMagn = _sqMagn; }
    /**
     * Return average vector of this vector and passed vector
     * @param q
     * @return
     */
    public myVectorf _avgWithMe(myVectorf q) {return new myVectorf((this.x+q.x)/2.0f,(this.y+q.y)/2.0f,(this.z+q.z)/2.0f);}
    /**
     * Static Method : Return average of the two passed vectors
     * @param p
     * @param q
     * @return
     */
    public static myVectorf _average(myVectorf p, myVectorf q) {return new myVectorf((p.x+q.x)/2.0f,(p.y+q.y)/2.0f,(p.z+q.z)/2.0f);}
    /**
     * Multiply this vector by n. Returns this vector
     */
    @Override
    public myVectorf _mult(float n){ super._mult(n); this._mag(); return this; }
    /**
     * Multiply this vector by n. Returns this vector
     */
    public myVectorf _mult(double n){ super._mult((float)n); this._mag(); return this; }
    /**
     * Static Method : Multiply passed vector p by n and returns result
     * @param p
     * @param n
     * @return
     */
    public static myVectorf _mult(myVectorf p, float n){ myVectorf result = new myVectorf(p.x * n, p.y * n, p.z * n); return result;}
    /**
     * Static Method : Multiply passed vector p by n and returns result
     * @param p
     * @param n
     * @return
     */
    public static myVectorf _mult(myVectorf p, double n){ myVectorf result = new myVectorf(p.x * n, p.y * n, p.z * n); return result;}
    /**
     * Static Method : Elementwise-multiply passed vector p by passed vector q and returns resultant vector
     * @param p
     * @param q
     * @return
     */
    public static myVectorf _mult(myVectorf p, myVectorf q){ myVectorf result = new myVectorf(p.x *q.x, p.y * q.y, p.z * q.z); return result;}    
    /**
     * Static Method : Elementwise-multiply passed vector p by passed vector q and returns resultant vector
     * @param p
     * @param q
     * @return
     */
    public static myVectorf _ewise_mult(myVectorf p, myVectorf q){ myVectorf result = new myVectorf(p.x *q.x, p.y * q.y, p.z * q.z); return result;}    
    /**
     * Static Method : Elementwise-multiply first two arguments and put result in 3rd argument
     * @param p
     * @param q
     * @param r
     */
    public static void _mult(myVectorf p, myVectorf q, myVectorf r){ myVectorf result = new myVectorf(p.x *q.x, p.y * q.y, p.z * q.z); r.set(result);}  
    /**
     * Divide this vector by q
     */
    @Override
    public void _div(float q){super._div(q); this._mag();}  
    /**
     * Static Method : Divide vector p by n
     * @param p
     * @param n
     * @return
     */
    public static myVectorf _div(myVectorf p, float n){ if(n==0) return p; myVectorf result = new myVectorf(p.x / n, p.y / n, p.z / n); return result;}
    /**
     * Static Method : Element-wise divide vector p by vector q
     * @param p
     * @param q
     * @return
     */
    public static myVectorf _div(myVectorf p, myVectorf q){ if((q.x==0)||(q.y==0)|| (q.z==0))return p; return new myVectorf(p.x / q.x, p.y / q.y, p.z / q.z); }
    /**
     * Static Method : Element-wise divide vector p by vector q
     * @param p
     * @param q
     * @return
     */
    public static myVectorf _ewise_div(myVectorf p, myVectorf q){ if((q.x==0)||(q.y==0)|| (q.z==0))return p; return new myVectorf(p.x / q.x, p.y / q.y, p.z / q.z); }
    /**
     * Add coordinate values to this vector
     */
    @Override
    public void _add(float _x, float _y, float _z){ super._add(_x, _y, _z); this._mag(); }
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
    public static myVectorf _add(myVectorf p, myVectorf q){ myVectorf result = new myVectorf(p.x + q.x, p.y + q.y, p.z + q.z); return result;}
    /**
     * Static Method : add vector to passed point and return result
     * @param p point
     * @param q vector
     * @return resulting point of element-wise addition of p + q
     */
    public static myVectorf _add(myPointf p, myVectorf q){ myVectorf result = new myVectorf(p.x + q.x, p.y + q.y, p.z + q.z); return result;}
    /**
     * Static Method : add two vectors, putting result in 3rd argument
     * @param p,q : vectors to add
     * @param r : resulting point of element-wise addition of p + q
     */
    public static void _add(myVectorf p, myVectorf q, myVectorf r){ myVectorf result = new myVectorf(p.x + q.x, p.y + q.y, p.z + q.z); r.set(result);}  
    /**
     * Static Method : Find the sum of an array of vectors, minimizes extra calcs.
     * @param vAra
     * @return
     */
    public static myVectorf _add(myVectorf[] vAra) {
        // aggregate to minimize magnitude calcs
        float x = 0, y = 0, z = 0;
        for(int i=0;i<vAra.length;++i) {
            x += vAra[i].x;
            y += vAra[i].y;
            z += vAra[i].z;
        }
        return new myVectorf(x,y,z);
    }
    /**
     * Subtract values from this vector
     */
    @Override
    public void _sub(float _x, float _y, float _z){ super._sub(_x, _y, _z);  this._mag(); }
    /**
     * Element-wise subtract passed vector from this vector
     * @param v
     */
    public void _sub(myVectorf v){ this.x -= v.x; this.y -= v.y; this.z -= v.z;  this._mag(); }
    /**
     * Static Method : Subtract q from p and return result
     * @param p
     * @param q
     * @return
     */
    public static myVectorf _sub(myPointf p, myPointf q){ return new myVectorf(p.x - q.x, p.y - q.y, p.z - q.z);}
    /**
     * Static Method : Subtract q from p, putting result in 3rd argument
     * @param p 
     * @param q
     * @param r
     */
    public static void _sub(myVectorf p, myVectorf q, myVectorf r){ myVectorf result = new myVectorf(p.x - q.x, p.y - q.y, p.z - q.z); r.set(result);}  
    /**
     * Calc, set and return this vector's magnitude
     * @return
     */
    public float _mag(){ this.magn = (float)Math.sqrt(this._SqMag()); return this.magn; }  
    /**
     * Calc, set and return this vector's sq magnitude
     * @return
     */
    public float _SqMag(){ this.sqMagn = ((this.x*this.x) + (this.y*this.y) + (this.z*this.z)); return this.sqMagn; }
    /**
     * Sets length of vector
     * @param _newMag
     */
    public void _scale(float _newMag){this._normalize()._mult(_newMag);}
    /**
     * Normalize this vector and return it
     * @return
     */
    public myVectorf _normalize(){this._mag();if(magn==0){return this;}this.x /= magn; this.y /= magn; this.z /= magn; _mag();return this;}
    /**
     * Static Method : build a normalizeed version of the passed vector
     * @param v
     * @return
     */
    public static myVectorf _normalize(myVectorf v){double magn = v._mag(); if(magn==0){return v;} myVectorf newVec = new myVectorf( v.x / magn, v.y / magn, v.z / magn);return newVec;}// newVec._mag(); return newVec;}
    /**
     * Recalculate and return this vector's magnitude
     * @return
     */
    public float _norm(){return _mag();}
    /**
     * Static Method : Calculate and return the L2 norm of the passed vector
     * @param v
     * @return
     */
    public static float _L2Norm(myVectorf v){return (float)Math.sqrt(v._SqMag());}
    /**
     * Static Method : Calculate and return the squared L2 norm of the passed vector
     * @param v
     * @return
     */
    public static float _L2SqNorm(myVectorf v){return v._SqMag();}
    /**
     * Build new normalized version of this vector (this vector remains unchanged)
     * @return
     */
    public myVectorf _normalized(){float magn = this._mag(); myVectorf newVec = (magn == 0) ? (new myVectorf(0,0,0)) : (new myVectorf( this.x /magn, this.y / magn, this.z / magn)); newVec._mag(); return newVec;}
    /**
     * Build a copy of this vector
     */
    @Override
    public myVectorf cloneMe(){myVectorf retVal = new myVectorf(this.x, this.y, this.z); return retVal;}  
    /**
     * Static Method : find the squared distance between the two passed vectors 
     * @param q
     * @param r
     * @return
     */
    public static float _SqrDist(myVectorf q, myVectorf r){  return ((r.x - q.x)*(r.x - q.x)) + ((r.y - q.y)*(r.y - q.y)) + ((r.z - q.z)*(r.z - q.z));}
    /**
     * Static Method : find the distance between the two passed vectors 
     * @param q
     * @param r
     * @return
     */
    public static float _dist(myVectorf q, myVectorf r){  return (float)Math.sqrt(((r.x - q.x) *(r.x - q.x)) + ((r.y - q.y) *(r.y - q.y)) + ((r.z - q.z) *(r.z - q.z)));}
    /**
     * Static Method : find the distance between the passed vector and the vector represented by the passed 3 values
     * @param r
     * @param qx
     * @param qy
     * @param qz
     * @return
     */
    public static float _dist(myVectorf r, float qx, float qy, float qz){  return (float)Math.sqrt(((r.x - qx) *(r.x - qx)) + ((r.y - qy) *(r.y - qy)) + ((r.z - qz) *(r.z - qz)));}    
    /**
     * Find the cross product of this vectoor and the passed vector
     * @param b
     * @return
     */
    public myVectorf _cross(myVectorf b){ return new myVectorf((this.y * b.z) - (this.z*b.y), (this.z * b.x) - (this.x*b.z), (this.x * b.y) - (this.y*b.x));}
    /**
     * Static Method : Find the cross product of the two passed vectors
     * @param a
     * @param b
     * @return
     */
    public static myVectorf _cross(myVectorf a, myVectorf b){        return a._cross(b);}
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
    public static myVectorf _cross(float ax, float ay, float az, float bx, float by, float bz){        return new myVectorf((ay*bz)-(az*by), (az*bx)-(ax*bz), (ax*by)-(ay*bx));}
    /**
     * Find the dot product of this vector dotted with the passed vector
     * @param b
     * @return
     */
    public float _dot(myVectorf b){return ((this.x * b.x) + (this.y * b.y) + (this.z * b.z));}
    /**
     * Static Method : Find the dot product of two two passed vectors
     * @param a
     * @param b
     * @return
     */
    public static float _dot(myVectorf a, myVectorf b){        return a._dot(b);}
    /**
     * Static Method : U|V det product
     * @param U
     * @param V
     * @return
     */
    public static float _det3(myVectorf U, myVectorf V) {float udv = U._dot(V); return (float)(Math.sqrt(U._dot(U)*V._dot(V) - (udv*udv))); }
    /**
     * Static Method : U*(VxW)  mixed product, determinant - measures 6x the volume of the parallelapiped formed by passed vectors
     * @param U
     * @param V
     * @param W
     * @return
     */
    public static float _mixProd(myVectorf U, myVectorf V, myVectorf W) {return U._dot(myVectorf._cross(V,W)); }
    /**
     * Static Method : U * (VxW)>0  U,V,W are clockwise
     * @param U
     * @param V
     * @param W
     * @return
     */
    public static boolean _isCW_Vecs(myVectorf U, myVectorf V, myVectorf W) {return _mixProd(U,V,W)>0; }
    /**
     * Static Method : Area of triangle described by 3 points 
     * @param A, B, C Triangle verts
     * @return area of proscribed triangle
     */
    public static float _area(myPointf A, myPointf B, myPointf C) {    myVectorf x = new myVectorf(A,B), y = new myVectorf(A,C), z = x._cross(y);     return z.magn/2.0f; } 

    /**
     * Static Method : Returns volume of tetrahedron defined by A,B,C,D
     * @param A,B,C,D verts of tet
     * @return volume
     */
    public static float _volume(myPointf A, myPointf B, myPointf C, myPointf D) {return _mixProd(new myVectorf(A,B),new myVectorf(A,C),new myVectorf(A,D))/6.0f; } 
    /**
     * Static Method : Returns true if tet is oriented so that A sees B, C, D clockwise
     * @param A,B,C,D verts of tet
     * @return if tet is oriented clockwise (A->B->C->D)
     */
    public static boolean _isCW_Tet(myPointf A, myPointf B, myPointf C, myPointf D) {return _volume(A,B,C,D)>0; } 
    /**
     * Static Method : Returns true if U and V are almost parallel
     * @param U
     * @param V
     * @return
     */
    public static boolean isParallel(myVectorf U, myVectorf V) {return U._cross(V).magn<U.magn*V.magn*MyMathUtils.EPS_F; }
    
    /**
     * Static Method : This will return a unit normal ortho to the plane spanned by U and V. If U and V are parallel, will provide a normal to U.
     * @param U 
     * @param V
     * @return
     */
    public static myVectorf _findNormToPlane(myVectorf U, myVectorf V) {
        myVectorf norm = V._cross(U);
        if(norm.magn<U.magn*V.magn*MyMathUtils.EPS_F) {
            //parallel, find normal to U and modified V 
            //that will never be coincident no matter what V is
            norm = myVectorf._cross(U.x, U.y, U.z, (V.x*2)+1, (V.y*3)+2, (V.z*4)+3);            
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
    public static float _angleBetween(myVectorf v1, myVectorf v2) {
        float  _v1Mag = v1._mag(), 
                _v2Mag = v2._mag(), 
                dotProd = v1._dot(v2),
                cosAngle = dotProd/(_v1Mag * _v2Mag),
                angle = (float)(Math.acos(cosAngle));
        return angle;
    }//_angleBetween
    /**
     * Returns the angle between this and the passed vector
     * @param v2
     * @return
     */   
    public float angleWithMe(myVectorf v2) {
        float  _v1Mag = _mag(), 
                _v2Mag = v2._mag(), 
                dotProd = _dot(v2),
                cosAngle = dotProd/(_v1Mag * _v2Mag),
                angle = (float) Math.acos(cosAngle);
        return angle;
    }//_angleBetween
    /**
     * Static Method : Rotate v1 around axis unit vector u, by Pi/2, around origin
     * @param v1
     * @param u
     * @return
     */
    public static myVectorf _rotAroundAxis(myVectorf v1, myVectorf u){return _rotAroundAxis(v1, u, MyMathUtils.HALF_PI_F);}
    /**
     * Static Method : Rotate v1 around axis unit vector u, by given angle thet, around origin
     * @param v1
     * @param u
     * @param thet
     * @return
     */
    public static myVectorf _rotAroundAxis(myVectorf v1, myVectorf u, float thet){        
        float cThet = (float)(Math.cos(thet)), sThet = (float)(Math.sin(thet)), oneMC = 1-cThet,
                ux2 = u.x*u.x, uy2 = u.y*u.y, uz2 = u.z*u.z,
                uxy = u.x * u.y, uxz = u.x * u.z, uyz = u.y*u.z,
                uzS = u.z*sThet, uyS = u.y*sThet, uxS = u.x*sThet,
                uxzC1 = uxz *oneMC, uxyC1 = uxy*oneMC, uyzC1 = uyz*oneMC;
        //build rot matrix in vector def
        myVectorf res = new myVectorf(
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
    public myVectorf rotMeAroundAxis(myVectorf u, double thet){        
        double cThet = Math.cos(thet), sThet = Math.sin(thet), oneMC = 1-cThet,
                ux2 = u.x*u.x, uy2 = u.y*u.y, uz2 = u.z*u.z,
                uxy = u.x * u.y, uxz = u.x * u.z, uyz = u.y*u.z,
                uzS = u.z*sThet, uyS = u.y*sThet, uxS = u.x*sThet,
                uxzC1 = uxz *oneMC, uxyC1 = uxy*oneMC, uyzC1 = uyz*oneMC;
        //build rot matrix in vector def
        myVectorf res = new myVectorf(
                (ux2*oneMC+cThet) * this.x + (uxyC1-uzS)         * this.y + (uxzC1+uyS) *this.z,
                (uxyC1+uzS)       * this.x + (uy2*oneMC+cThet)* this.y + (uyzC1-uxS) *this.z,
                (uxzC1-uyS)       * this.x + (uyzC1+uxS)        * this.y + (uz2*oneMC + cThet) * this.z);
        
        return res;        
    }

    /**
     * undefined for a vector - use rotMeAroundAxis
     */
    public final myPointf rotMeAroundPt(float a, myVectorf I, myVectorf J, myPointf G) {return this;}
    /**
     * undefined for a vector - rotMeAroundAxis
     */
    public final myPointf rotMeAroundPt(myPointf C, myPointf P, myPointf R)  {return this;}
    /**
     * Static Method :  Find the angle between two vectors - Note this version is for 2D - relies on neither vector being coplanar with (0,0,1);
     * @param U
     * @param V
     * @return
     */    
    public static float _angleBetween_Xprod(myVectorf U, myVectorf V){
        return _angleBetween_Xprod(U,V, new myVectorf(0,0,1));
    }

    /**
     * Static Method : Find the angle between two vectors, with axis about which to determine sign
     * @param U
     * @param V
     * @param axis unit length axis about which to determine sign
     * @return
     */    
    public static float _angleBetween_Xprod(myVectorf U, myVectorf V, myVectorf axis){
        myVectorf cross = U._cross(V);
        double dot = U._dot(V);
        
        float angle = (float) Math.atan2(cross.magn,dot),
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
        if (b instanceof myVectorf v) {return ((this.x == v.x) && (this.y == v.y) && (this.z == v.z));}
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
}//myVectorf




