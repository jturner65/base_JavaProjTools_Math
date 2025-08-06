package base_Math_Objects.vectorObjs.floats;

public class myQuaternionf {
    /**
     * Vector representation/direction
     */
    myVectorf v = new myVectorf();
    /**
     * Scalar quantities
     */
    public float w,  magn, sqMagn;
    /**
     * Empty unit quat
     */
    public myQuaternionf(){w = 1.0f;_mag();}
    /**
     * Vector + scalar representation
     * @param _v
     * @param _w
     */
    public myQuaternionf(myVectorf _v, float _w){v = new myVectorf(_v);  w = _w; _mag();    } 
    /**
     * Builds vector representation from _x,_y,_z
     * @param _x
     * @param _y
     * @param _z
     * @param _w
     */
    public myQuaternionf(float _x, float _y, float _z, float _w) {v = new myVectorf(_x,_y,_z); w = _w; _mag();    }
    /**
     * Builds vector representation from _x,_y,_z
     * @param _x
     * @param _y
     * @param _z
     * @param _w
     */
    public myQuaternionf(double _x, double _y, double _z, double _w) {v = new myVectorf(_x,_y,_z); w = (float) _w; _mag();   }
    /**
     * Copy constructor
     * @param a
     */
    public myQuaternionf(myQuaternionf a) {this(new myVectorf(a.v),a.w);    }    
    
    //call internally whenever values change to keep vector and x,y,z values synched
    protected void _pset(float _x, float _y, float _z, float _w){v.set(_x,_y,_z); w = _w;_mag();}
    protected void _pset(myVectorf _v, float _w){v.set(_v); w = _w;_mag();}    
    /**
     * given axis angle representation, convert to quaternion
     * @param theta
     * @param vec
     */
    public void setFromAxisAngle(float theta, myVectorf vec){
        vec._normalize();
        float htht = theta/2.0f, sThetH = (float)(Math.sin(htht));
        _pset(myVectorf._mult(vec,sThetH), (float)(Math.cos(htht)));
    }
    //
    public float _mag(){ this.magn = (float)Math.sqrt(this._SqMag()); return magn; }  
    public float _SqMag(){ this.sqMagn = this.v.sqMagn + (this.w*this.w); return this.sqMagn; }                             //squared magnitude
    public float _dot(myQuaternionf q){return (this.v._dot(q.v) + (this.w*q.w));}

    public myQuaternionf _normalize(){this._mag();if(magn==0){return this;}_pset(v.x/magn,v.y/magn,v.z/magn,w/magn);return this;}
    public static myQuaternionf _normalize(myQuaternionf _q){_q._mag(); return new myQuaternionf(_q.v.x/_q.magn,_q.v.y/_q.magn,_q.v.z/_q.magn,_q.w/_q.magn);}
    //multiply this quaternion by another quaternion -> this * q
    public myQuaternionf _qmult(myQuaternionf q){
        //w1v2 + w2v1 + v1.cross(v2)
        myVectorf t1 = myVectorf._mult(q.v, w), t2 = myVectorf._mult(v, q.w), t3 = v._cross(q.v);
        t1._add(t2);t1._add(t3);
        myQuaternionf res = new myQuaternionf(t1, (this.w * q.w) - v._dot(q.v) );
        return res;
    }
    //give the conjugate of this quaternion
    public myQuaternionf _conj(){return new myQuaternionf(myVectorf._mult(v, -1.0f),w);}
    
    //rotate the passed vector by this quaternion -> q*q_v*qstar
    public myVectorf _rot(myVectorf _v){        
        myQuaternionf qnorm = myQuaternionf._normalize(this),
                q_v = new myQuaternionf(_v, 0),
                res1 = qnorm._qmult(q_v),
                conj = qnorm._conj(), 
                res = res1._qmult(conj);
        return res.v;
    }
    
    //rotate toRotVec by passed angle around rVec
    public static myVectorf _quatRot(float theta, myVectorf rVec, myVectorf toRotVec){
        myQuaternionf _tmp = new myQuaternionf();
        _tmp.setFromAxisAngle(theta,rVec);
        return _tmp._rot(toRotVec);    
    }
    
    //convert this to axis angle - theta,rx,ry,rz
    public float[] toAxisAngle(){
        myQuaternionf tmp = new myQuaternionf(this);
        if (tmp.w > 1) {tmp._normalize();}                                                         
        float thet = 2 * (float)Math.acos(tmp.w),s = (float)Math.sqrt(1-tmp.w*tmp.w);             
        if (s < 0.0000001) { return new float[]{thet, 1,0,0};    }            //with s close to 0, thet doesn't matter 
        else {            return new float[]{thet, tmp.v.x/s, tmp.v.y/s, tmp.v.z/s};   }
    }//asAxisAngle
    
    private static myQuaternionf _lerp(myQuaternionf qa, myQuaternionf qb, float t1, float t2){return new myQuaternionf(myVectorf._add(myVectorf._mult(qa.v,t1), myVectorf._mult(qb.v,t2)), (qa.w * t1) + (qb.w*t2));}
    
    public static myQuaternionf _slerp(myQuaternionf qa, myQuaternionf qb, float t){
        // Calculate angle between them.
        float cosHalfTheta = qa._dot(qb);
        // if qa=qb or qa=-qb then theta = 0 and we can return qa
        if (cosHalfTheta >= 1.0){    return new myQuaternionf(qa);    }
        float halfTheta = (float) Math.acos(cosHalfTheta);
        float s = (float) Math.sqrt(1.0 - cosHalfTheta*cosHalfTheta);
        if (s < 0.0000001){    return _lerp(qa, qb,.5f, .5f);}//avoid div by zero - any t will do since denotes axis
        return _lerp(qa, qb, (float) Math.sin((1 - t) * halfTheta)/s , (float) Math.sin(t * halfTheta) / s);
    }//_slerp
    
    /**
     * returns if this quaternion is equal to passed quaternion
     * @param b myQuaternion to check
     * @return whether they are equal
     */
    @Override
    public boolean equals(Object b){
        if (this == b) return true;
        if (b instanceof myQuaternionf q) {return ((this.v.x == q.v.x) && (this.v.y == q.v.y) && (this.v.z == q.v.z) && (this.w == q.w));}
        return false;
    }
    
    public String toString(){return "vec:"+v.toStrBrf() + "\t w:"+ String.format("%.4f",w);} 
    
}//myQuaternionf