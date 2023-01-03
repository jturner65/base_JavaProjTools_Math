package base_Math_Objects;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;

import base_Math_Objects.vectorObjs.doubles.myPoint;
import base_Math_Objects.vectorObjs.floats.myPointf;
import base_Math_Objects.vectorObjs.doubles.myVector;
import base_Math_Objects.vectorObjs.floats.myVectorf;

/**
 * mathematical functions and constants that might be of use in some applications
 * @author john
 *
 */
public class MyMathUtils {
	/**
	 * Useful double constants
	 */
	public static final double 
		PI = Math.PI,
		QUARTER_PI = .25*PI,
		HALF_PI = .5*PI,
		THIRD_PI = PI/3.0,
		TWO_PI = 2.0*PI,
		THREE_QTR_PI = .75 * PI,
		FIFTH_PI = .2 * PI,
		SQRT_2 = Math.sqrt(2.0),
		INV_SQRT_2 = .5 * SQRT_2,
		SQRT_3 = Math.sqrt(3.0),
		INV_SQRT_3 = 1.0/SQRT_3,
		DEG_TO_RAD = PI/180.0,
		EPS = 1e-8,
		LOG_2 = Math.log(2.0),
		LOG_10 = Math.log(10.0);
	
	/**
	 * Useful float constants
	 */
	public static final float 
		PI_F = (float)PI,
		QUARTER_PI_F = (float)QUARTER_PI,
		HALF_PI_F = (float) HALF_PI,
		TWO_PI_F = (float) TWO_PI,
		THIRD_PI_F = (float) THIRD_PI,
		THREE_QTR_PI_F = (float) THREE_QTR_PI,
		FIFTH_PI_F = (float) FIFTH_PI,
		SQRT_2_F = (float) SQRT_2,
		INV_SQRT_2_F = (float) INV_SQRT_2,
		SQRT_3_F = (float) SQRT_3,
		INV_SQRT_3_F = (float) INV_SQRT_3,
		DEG_TO_RAD_F = (float)DEG_TO_RAD,
		EPS_F = (float) EPS,
		LOG_2_F = (float)LOG_2,
		LOG_10_F = (float)LOG_10;
	

    // numbers greater than 10^MAX_DIGITS_10 or e^MAX_DIGITS_EXP are considered unsafe ('too big') for floating point operations
    protected static final int MAX_DIGITS_EXP = 677;
    protected static final int MAX_DIGITS_10 = 294; // ~ MAX_DIGITS_EXP/LN(10)
    protected static final int MAX_DIGITS_2 = 977; // ~ MAX_DIGITS_EXP/LN(2)
        
    /**
     * # of pre-computed cosine and sine values.  
     */
    public static final float numTrigVals = 36.0f;
    /**
     * Pre-computed array of {@value #numTrigVals} cosine values 0->2pi
     */
    public static final double[] preCalcCosVals;
    /**
     * Pre-computed array of {@value #numTrigVals} sine values 0->2pi
     */
    public static final double[] preCalcSinVals;
    
    private static final float deltaThet = TWO_PI_F/numTrigVals,
			finalThet = TWO_PI_F+deltaThet;
    /**
     * Pre-computed array of {@value #numTrigVals} float cosine values 0->2pi
     */
    public static final float[] preCalcCosVals_f;
    /**
     * Pre-computed array of {@value #numTrigVals} float sine values 0->2pi
     */
    public static final float[] preCalcSinVals_f;
    
    /**
     * Build arrays of precomputed sine and cosine values
     */
    static {
      	preCalcCosVals = new double[(int)numTrigVals + 2];
      	preCalcSinVals = new double[(int)numTrigVals + 2];
	   	preCalcCosVals_f = new float[(int)numTrigVals + 2];
	   	preCalcSinVals_f = new float[(int)numTrigVals + 2];
		int i=0;
		for(float a=0; a<=finalThet; a+=deltaThet) {
			preCalcCosVals[i] = Math.cos(a);
			preCalcSinVals[i] = Math.sin(a);
			preCalcCosVals_f[i] = (float) preCalcCosVals[i];
			preCalcSinVals_f[i] = (float) preCalcSinVals[i];
			++i;
		}      	
     }

	//shouldn't be instanced
	private MyMathUtils() {	}
	
	/**
	 * find distance of P to line determined by AB
	 * @param P
	 * @param A
	 * @param B
	 * @return
	 */
	public static double distToLine(myPoint P, myPoint A, myPoint B) {
		myVector AB = new myVector(A,B),AP = new myVector(A,P);
		AB._normalize();
		double udv = AB._dot(AP); //project AP onto line		
		double resSq = AP.sqMagn - (udv*udv);		//AB mag is 1 so can be ignored
		if(resSq <= 0) {return 0;}
		return Math.sqrt(resSq); 
	}
	
	/**
	 * find distance of P to line determined by AB
	 * @param P
	 * @param A
	 * @param B
	 * @return
	 */	
	public static float distToLine(myPointf P, myPointf A, myPointf B) {
		myVectorf AB = new myVectorf(A,B),AP = new myVectorf(A,P);
		AB._normalize();
		double udv = AB._dot(AP); //project AP onto line		
		double resSq = AP.sqMagn - (udv*udv);		//AB mag is 1 so can be ignored
		if(resSq <= 0) {return 0;}
		return (float) Math.sqrt(resSq); 
	}
	
	/**
	 * return the projection point of P on line determined by AB between A and B
	 * @param P point to investigate
	 * @param A,B line seg endpoints
	 * @return 
	 */	
	public static myPoint projectionOnLine(myPoint P, myPoint A, myPoint B) {
		myVector AB = new myVector(A,B), AP = new myVector(A,P);
		return new myPoint(A,AB._dot(AP)/(AB._dot(AB)),AB);
	}
	/**
	 * return the projection point of P on line determined by AB between A and B
	 * @param P point to investigate
	 * @param A,B line seg endpoints
	 * @return 
	 */	
	public static myPointf projectionOnLine(myPointf P, myPointf A, myPointf B) {
		myVectorf AB = new myVectorf(A,B), AP = new myVectorf(A,P);
		return new myPointf(A,AB._dot(AP)/(AB._dot(AB)),AB);
	}
	/**
	 * return true if P orthogonally projects onto line determined by AB between A and B
	 * @param P point to investigate
	 * @param A,B line seg endpoints
	 * @return
	 */
	public static boolean projectsBetween(myPoint P, myPoint A, myPoint B) {
		myVector AP = new myVector(A,P), AB = new myVector(A,B);
		if(AP._dot(AB) <= 0) {return false;}		//if not greater than 0 than won't project onto AB - past A away from segment
		myVector BP = new myVector(B,P), BA = new myVector(B,A);		
		return BP._dot(BA)>0 ; 						//if not greater than 0 than won't project onto AB - past B away from segment
	}
	/**
	 * return true if P orthogonally projects onto line determined by AB between A and B
	 * @param P point to investigate
	 * @param A,B line seg endpoints
	 * @return
	 */
	public static boolean projectsBetween(myPointf P, myPointf A, myPointf B) {
		myVectorf AP = new myVectorf(A,P), AB = new myVectorf(A,B);
		if(AP._dot(AB) <= 0) {return false;}		//if not greater than 0 than won't project onto AB - past A away from segment
		myVectorf BP = new myVectorf(B,P), BA = new myVectorf(B,A);		
		return BP._dot(BA)>0 ; 						//if not greater than 0 than won't project onto AB - past B away from segment
	}
	
	//public static final int O_FWD = 0, O_RHT = 1,  O_UP = 2;
	/**
	 * build axis angle orientation from passed orientation matrix
	 * @param orientation array of 3 vectors corresponding to orientation vectors
	 * @param O_FWD idx of forward orientation
	 * @param O_RHT idx of right orientation
	 * @param O_UP idx of up orientation
	 * @return axis-angle representation of orientation
	 */
	public static float[] toAxisAngle(myVectorf[] orientation, int O_FWD, int O_RHT, int O_UP) {
		float rt2 = .5f*SQRT_2_F;
		float angle,x=rt2,y=rt2,z=rt2,s;
		float fyrx = -orientation[O_FWD].y+orientation[O_RHT].x,
			uxfz = -orientation[O_UP].x+orientation[O_FWD].z,
			rzuy = -orientation[O_RHT].z+orientation[O_UP].y;
		float epsValCalcSq = EPS_F*EPS_F;
		if (((fyrx*fyrx) < epsValCalcSq) && ((uxfz*uxfz) < epsValCalcSq) && ((rzuy*rzuy) < epsValCalcSq)) {			//checking for rotational singularity
			// angle == 0
			float fyrx2 = orientation[O_FWD].y+orientation[O_RHT].x,
				fzux2 = orientation[O_FWD].z+orientation[O_UP].x,
				rzuy2 = orientation[O_RHT].z+orientation[O_UP].y,
				fxryuz3 = orientation[O_FWD].x+orientation[O_RHT].y+orientation[O_UP].z-3;
			if (((fyrx2*fyrx2) < 1)	&& (fzux2*fzux2 < 1) && ((rzuy2*rzuy2) < 1) && ((fxryuz3*fxryuz3) < 1)) {	return new float[]{0,1,0,0}; }
			// angle == pi
			angle = PI_F;
			float fwd2x = (orientation[O_FWD].x+1)/2.0f,rht2y = (orientation[O_RHT].y+1)/2.0f,up2z = (orientation[O_UP].z+1)/2.0f,
				fwd2y = fyrx2/4.0f, fwd2z = fzux2/4.0f, rht2z = rzuy2/4.0f;
			if ((fwd2x > rht2y) && (fwd2x > up2z)) { // orientation[O_FWD].x is the largest diagonal term
				if (fwd2x< EPS_F) {	x = 0;} else {			x = (float) Math.sqrt(fwd2x);y = fwd2y/x;z = fwd2z/x;} 
			} else if (rht2y > up2z) { 		// orientation[O_RHT].y is the largest diagonal term
				if (rht2y< EPS_F) {	y = 0;} else {			y = (float) Math.sqrt(rht2y);x = fwd2y/y;z = rht2z/y;}
			} else { // orientation[O_UP].z is the largest diagonal term so base result on this
				if (up2z< EPS_F) {	z = 0;} else {			z = (float) Math.sqrt(up2z);	x = fwd2z/z;y = rht2z/z;}
			}
			return new float[]{angle,x,y,z}; // return 180 deg rotation
		}
		//no singularities - handle normally
		myVectorf tmp = new myVectorf(rzuy, uxfz, fyrx);
		s = tmp.magn;
		if (s < EPS_F){ s=1; }
		tmp._scale(s);//changes mag to s
			// prevent divide by zero, should not happen if matrix is orthogonal -- should be caught by singularity test above
		angle = (float) -Math.acos(( orientation[O_FWD].x + orientation[O_RHT].y + orientation[O_UP].z - 1)/2.0);
		
		//consume this as follows : 
		//p.rotate(O_axisAngle[0],O_axisAngle[1],O_axisAngle[2],O_axisAngle[3]);
		
	   return new float[]{angle,tmp.x,tmp.y,tmp.z};
	}//toAxisAngle
	
	
	
	/**
	 * Calculate normal to planed described by triangle ABC, non-normalized (proportional to area)
	 * @param A, B, C verts of triangle
	 * @return
	 */
	public myVector normToPlane(myPoint A, myPoint B, myPoint C) {
		return myVector._cross(new myVector(A,B),new myVector(A,C)); 
	};   // normal to triangle (A,B,C), not normalized (proportional to area)
	
	/**
	 * Calculate normal to planed described by triangle ABC, non-normalized (proportional to area)
	 * @param A, B, C verts of triangle
	 * @return
	 */
	public myVectorf normToPlane(myPointf A, myPointf B, myPointf C) {
		return myVectorf._cross(new myVectorf(A,B),new myVectorf(A,C)); 
	};   // normal to triangle (A,B,C), not normalized (proportional to area)
	
	/**
	 * Kahan (quake) inv sqrt calc - can be up to 30% faster than Math lib
	 * @param x
	 * @return
	 */
    public static float invSqrtFloat(float x){
        float xhalf = x * 0.5f;
        int i = Float.floatToIntBits(x);
        //hex number is floating point/bitwise rep of approx sqrt(2^127)
        i = 0x5f3759df - (i >> 1);
        x = Float.intBitsToFloat(i);
        x *= (1.5f - (xhalf * x * x));   // newton iter
        return x;
    }
    /**
     * double version of Kahan (quake) sqrt approx - can be up to 30% faster than Math lib
     * @param x
     * @return
     */
    public static double invSqrtDouble(double x){    	
        double xhalf = x * 0.5;
        long i = Double.doubleToLongBits(x);
        //hex number is double/bitwise rep of approx sqrt(2^127)
        i = 0x5fe6eb50c7b537a9L - (i >> 1);
        x = Double.longBitsToDouble(i);
        x *= (1.5 - (xhalf * x * x));   // newton iter
        return x;
    }
    
    /**
     * iterative factorial formulation
     * @param x
     * @return
     */
    public static long fact(int x) {
    	long ttl=1;
    	for(int i=x; i>1;--i) { ttl*=i;}
    	return ttl;
    }
     
	/**
	 * Use column math to calculate x * values in ans ara
	 * @param x
	 * @param ans
	 * @param MSigDig
	 * @return
	 */
	private static int calcMult(int x, byte[] ans, int MSigDig) {
		long carry = 0, prod;
		
		for(int i=0; i<MSigDig; ++i) {			
			prod = ans[i] * x + carry;
			ans[i] = (byte) (prod % 10);
			carry = prod/10;
		}
		//propagate carry
	    while (carry>0) { 
	        ans[MSigDig] = (byte) (carry%10); 
	        carry = carry/10; 
	        MSigDig++; 
	    } 
	    return MSigDig; 
	}//calcMult
    
    /**
     * Find factorial of huge number, putting result in an array.
     * @param x
     * @return x! (factorial of x) as array of digits base 10
     */
    public static byte[] bigFact(int x) {
    	if (x < 2) {return new byte[] {1};}
    	//find # of digits by finding 1 + ceil(log10(x!))  
    	//== log10(x) + log10(x-1) ... + log10(2)
    	double sum = 0.0;
    	for(int i=2;i<=x;++i) {sum += Math.log10(i); 	}
    	int numDigits = (int) (Math.ceil(sum)); 
    	byte[] res = new byte[numDigits];
    	res[0] = 1;
    	int MSigDig = 1; 
    	for (int i=2; i<=x;++i) {
    		MSigDig = calcMult(i, res, MSigDig);
    	}   
    	//now reverse array
    	int lastIdx = numDigits-1;
    	for(int i=lastIdx; i>=numDigits/2.0f; --i) {
    		byte tmp = res[i];
    		res[i] = res[lastIdx-i];
    		res[lastIdx-i] = tmp;
    	}
    	return res;
    }//bigFact
    
    /**
     * n choose k == n!/(k! * (n-k)!) - ways to choose k items from a set of size n
     * @param n top - size of pop
     * @param k bottom - size of choice
     * @return
     */
    public static long choose(long n, long k) {		//entry point
    	if(k>n) {return 0;} if (k==n) {return 1;}
    	if(n-k < k) {		return  _choose(n,n-k);   	}
    	return _choose(n,k);
    }
    private static long _choose(long n, long k) {		//multiplicative formulation
    	long res = 1;
    	for(int i=1;i<=k;++i) {    		res *= (n+1-i);res /=i;    	}
    	return res;
    }

    /**
     * n choose k == n!/(k! * (n-k)!) - ways to choose k items from a set of size n
     * @param n top - size of pop
     * @param k bottom - size of choice
     * @return
     */
    public static BigInteger choose_BigInt(long n, long k) {		//entry point
    	if(k>n) {return BigInteger.ZERO;} if (k==n) {return BigInteger.ONE;}
    	if(n-k < k) {		return  _choose_BigInt(n,n-k);   	}
    	return _choose_BigInt(n,k);
    }
    private static BigInteger _choose_BigInt(long n, long k) {		//multiplicative formulation
    	BigInteger res = BigInteger.ONE;
    	for(int i=1;i<=k;++i) {    
    		res = res.multiply(BigInteger.valueOf(n+1-i));
    		res = res.divide(BigInteger.valueOf(i));
    	}
    	return res;
    }

    /**
     * Computes the natural logarithm of a BigInteger. 
     * Works for really big integers (practically unlimited), even when the argument 
     * falls outside the double range
     * Returns Nan if argument is negative, NEGATIVE_INFINITY if zero.
     * @param val Argument
     * @return Natural logarithm, as in Math.log()
     */
    public static double logBigInteger(BigInteger val) {
        if (val.signum() < 1) { return val.signum() < 0 ? Double.NaN : Double.NEGATIVE_INFINITY;}
        int blex = val.bitLength() - MAX_DIGITS_2; // any value in 60..1023 works ok here
        if (blex > 0) {
            val = val.shiftRight(blex);
            double res = Math.log(val.doubleValue());
            return res + blex * LOG_2;
        } else {        	return Math.log(val.doubleValue());        }
    }

    /**
     * Computes the natural logarithm of a BigDecimal. 
     * Works for really big (or really small) arguments, even outside the double range.
     * Returns Nan if argument is negative, NEGATIVE_INFINITY if zero.
    *
     * @param val Argument
     * @return Natural logarithm, as in <tt>Math.log()</tt>
     */
    public static double logBigDecimal(BigDecimal val) {
        if (val.signum() < 1) { return val.signum() < 0 ? Double.NaN : Double.NEGATIVE_INFINITY;}
        int digits = val.precision() - val.scale(); 
        if (digits < MAX_DIGITS_10 && digits > -MAX_DIGITS_10) {return Math.log(val.doubleValue());}
        else {            return logBigInteger(val.unscaledValue()) - val.scale() * LOG_10;}
    }

    /**
     * Computes the exponential function, returning a BigDecimal (precision ~ 16).       
     * Works for very big and very small exponents, even when the result 
     * falls outside the double range
     *
     * @param exponent Any finite value (infinite or Nan throws IllegalArgumentException)
     * @return The value of e (base of the natural logarithms) raised to the given exponent, as in Math.exp()
     */
    public static BigDecimal expBig(double exponent) {
        if (!Double.isFinite(exponent)) {throw new IllegalArgumentException("Infinite not accepted: " + exponent);}
        // e^b = e^(b2+c) = e^b2 2^t with e^c = 2^t 
        double bc = MAX_DIGITS_EXP;
        if (exponent < bc && exponent > -bc) {return new BigDecimal(Math.exp(exponent), MathContext.DECIMAL64);}
        boolean neg = false;
        if (exponent < 0) {            neg = true;            exponent = -exponent;        }
        double b2 = bc, c = exponent - bc;
        int t = (int) Math.ceil(c / LOG_10);
        c = t * LOG_10;
        b2 = exponent - c;
        if (neg) {          b2 = -b2;         t = -t;   }
        return new BigDecimal(Math.exp(b2), MathContext.DECIMAL64).movePointRight(t);
    }

    /**
     * Same as Math.pow(a,b) but returns a BigDecimal (precision ~ 16). 
     * Works even for outputs that fall outside the double range
     * 
     * The only limit is that b * log(a) does not overflow the double range 
     * 
     * @param a Base. Should be non-negative 
     * @param b Exponent. Should be finite (and non-negative if base is zero)
     * @return Returns the value of the first argument raised to the power of the second argument.
     */
    public static BigDecimal powBig(double a, double b) {
        if (!(Double.isFinite(a) && Double.isFinite(b)))
            throw new IllegalArgumentException(Double.isFinite(b) ? "base not finite: a=" + a : "exponent not finite: b=" + b);
        if (b == 0) {  return BigDecimal.ONE;}
        if (b == 1) {  return BigDecimal.valueOf(a);}
        if (a == 0) {
            if (b >= 0) { return BigDecimal.ZERO;}
            else {        throw new IllegalArgumentException("0**negative = infinite");}
        }
        if (a < 0) { throw new IllegalArgumentException("negative base a=" + a);}
        double x = b * Math.log(a);
        if (Math.abs(x) < MAX_DIGITS_EXP) { return BigDecimal.valueOf(Math.pow(a, b));}
        else {          				    return expBig(x);}
    }
    
	/**
	 * calculate the normal, tangent, binormal components of passed vector compared to the passed normal (needs to be normalized)
	 * @param vec
	 * @param norm
	 * @return
	 */
	public static myVectorf[] getVecFrameNonNorm(myVectorf vec, myVectorf norm) {
		myVectorf[] result = new myVectorf[3];//(2, myVector(0, 0, 0));
		result[0] = myVectorf._mult(norm,(norm._dot(vec)));//norm dir
		result[1] = myVectorf._sub(vec, result[0]);		//tan dir
		result[2] = myVectorf._cross(result[0], result[1]);
		return result;
	}
    
	/**
	 * calculate the normal, tangent, binormal components of passed vector compared to the passed normal
	 * @param vec
	 * @param norm
	 * @return
	 */
	public static myVectorf[] getVecFrameNormalized(myVectorf vec, myVectorf norm) {
		myVectorf[] nn_result = getVecFrameNonNorm(vec, norm), result = new myVectorf[nn_result.length];
		for(int i=0;i<result.length;++i) {
			result[i]=nn_result[i]._normalized();
		}
		return result;
	}
    
	/**
	 * calculate the normal, tangent, binormal components of passed vector compared to the passed normal (needs to be normalized)
	 * @param vec
	 * @param norm
	 * @return
	 */
	public static myVector[] getVecFrameNonNorm(myVector vec, myVector norm) {
		myVector[] result = new myVector[3];//(2, myVector(0, 0, 0));
		result[0] = myVector._mult(norm,(norm._dot(vec)));//norm dir
		result[1] = myVector._sub(vec, result[0]);		//tan dir
		result[2] = myVector._cross(result[0], result[1]);
		return result;
	}
    
	/**
	 * calculate the normal, tangent, binormal components of passed vector compared to the passed normal
	 * @param vec
	 * @param norm
	 * @return
	 */
	public static myVector[] getVecFrameNormalized(myVector vec, myVector norm) {
		myVector[] nn_result = getVecFrameNonNorm(vec, norm), result = new myVector[nn_result.length];
		for(int i=0;i<result.length;++i) {
			result[i]=nn_result[i]._normalized();
		}
		return result;
	}
	
	/**
	 * Floor of passed float
	 * @param x
	 * @return
	 */
	public static int floor(float x){	return x>0 ? (int)x : (int)x-1;} 
	 
	/**
	 * Floor of passed double
	 * @param x
	 * @return
	 */
	public static int floor(double x){ return x>0 ? (int)x : (int)x-1;} 
	
	/**
	 * Return max of 2 ints
	 * @param x
	 * @param y
	 * @return
	 */
	public static int max(int x, int y) {      return (x>y) ? x : y;    }
	/**
	 * Return min of 2 ints
	 * @param x
	 * @param y
	 * @return
	 */
	public static int min(int x, int y) {      return (x<y) ? x : y;    }
	
	/**
	 * Return max of 2 longs
	 * @param x
	 * @param y
	 * @return
	 */
	public static long max(long x, long y) {      return (x>y) ? x : y;    }
	/**
	 * Return min of 2 longs
	 * @param x
	 * @param y
	 * @return
	 */
	public static long min(long x, long y) {      return (x<y) ? x : y;    }
	
	
	/**
	 * Return max of 2 floats
	 * @param x
	 * @param y
	 * @return
	 */
	public static float max(float x, float y) {      return (x>y) ? x : y;    }
	/**
	 * Return min of 2 floats
	 * @param x
	 * @param y
	 * @return
	 */
	public static float min(float x, float y) {      return (x<y) ? x : y;    }
	
	/**
	 * Return max of 2 doubles
	 * @param x
	 * @param y
	 * @return
	 */
	public static double max(double x, double y) {      return (x>y) ? x : y;    }
	/**
	 * Return min of 2 doubles
	 * @param x
	 * @param y
	 * @return
	 */
	public static double min(double x, double y) {      return (x<y) ? x : y;    }
	
	/**
	 * Return max of 3 ints
	 */
	public static int max(int x, int y, int z) {      return max(max(x,y),z);  }
	/**
	 * Return min of 3 ints
	 */
	public static int min(int x, int y, int z) {      return min(min(x,y),z);  }
	/**
	 * Return max of 3 longs
	 */
	public static long max(long x, long y, long z) {      return max(max(x,y),z);  }
	/**
	 * Return min of 3 longs
	 */
	public static long min(long x, long y, long z) {      return min(min(x,y),z);  }
	/**
	 * Return max of 3 floats
	 */
	public static float max(float x, float y, float z) {      return max(max(x,y),z);  }
	/**
	 * Return min of 3 floats
	 */
	public static float min(float x, float y, float z) {      return min(min(x,y),z);  }	
	/**
	 * Return max of 3 doubles
	 */
	public static double max(double x, double y, double z) {      return max(max(x,y),z);  }
	/**
	 * Return min of 3 doubles
	 */
	public static double min(double x, double y, double z) {      return min(min(x,y),z);  }
	
	
	/**
	 * Return max of array of ints
	 */
	public static int max(int[] valAra) {      
		int res = valAra[0];
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] > res){res = valAra[i];}	
		}
		return res; 
	}
	/**
	 * Return min of array of ints
	 */
	public static int min(int[] valAra) {
		int res = valAra[0];
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] < res){res = valAra[i];}	
		}
		return res; 
	}
	/**
	 * Return max of array of longs
	 */
	public static long max(long[] valAra) {      
		long res = valAra[0];
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] > res){res = valAra[i];}	
		}
		return res; 
	}
	/**
	 * Return min of array of longs
	 */
	public static long min(long[] valAra) {
		long res = valAra[0];
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] < res){res = valAra[i];}	
		}
		return res; 
	}
	/**
	 * Return max of array of floats
	 */
	public static float max(float[] valAra) {
		float res = valAra[0];
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] > res){res = valAra[i];}	
		}
		return res; 
	}
	/**
	 * Return min of array of floats
	 */
	public static float min(float[] valAra) {
		float res = valAra[0];
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] < res){res = valAra[i];}	
		}
		return res; 
	}	
	/**
	 * Return max of array of double
	 */
	public static double max(double[] valAra) {
		double res = valAra[0];
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] > res){res = valAra[i];}	
		}
		return res; 
	}
	/**
	 * Return min of array of double
	 */
	public static double min(double[] valAra) {
		double res = valAra[0];
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] < res){res = valAra[i];}	
		}
		return res; 
	}
	
	/**
	 * Return min(idx 0) and max (idx 1) of passed array of int values
	 * @param valAra one or more values 
	 * @return
	 */
	public static int[] minAndMax(int[] valAra) {    		
		int[] res = new int[] {valAra[0], valAra[0]};
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] < res[0]){res[0] = valAra[i];}	//min value
			if(valAra[i] > res[1]){res[1] = valAra[i];}	//max value
		}
		return res;
	}
	
	/**
	 * Return min(idx 0) and max (idx 1) of passed array of long values
	 * @param valAra one or more values 
	 * @return
	 */
	public static long[] minAndMax(long[] valAra) {    		
		long[] res = new long[] {valAra[0], valAra[0]};
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] < res[0]){res[0] = valAra[i];}	//min value
			if(valAra[i] > res[1]){res[1] = valAra[i];}	//max value
		}
		return res;
	}
	
	/**
	 * Return min(idx 0) and max (idx 1) of passed array of float values
	 * @param valAra one or more values 
	 * @return
	 */
	public static float[] minAndMax(float[] valAra) {    		
		float[] res = new float[] {valAra[0], valAra[0]};
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] < res[0]){res[0] = valAra[i];}	//min value
			if(valAra[i] > res[1]){res[1] = valAra[i];}	//max value
		}
		return res;
	}
	
	/**
	 * Return min(idx 0) and max (idx 1) of passed array of double values
	 * @param valAra one or more values 
	 * @return
	 */
	public static double[] minAndMax(double[] valAra) {    		
		double[] res = new double[] {valAra[0], valAra[0]};
		for (int i=1;i<valAra.length;++i) {	
			if(valAra[i] < res[0]){res[0] = valAra[i];}	//min value
			if(valAra[i] > res[1]){res[1] = valAra[i];}	//max value
		}
		return res;
	}
	
	
	/**
     * return max value of any comparable type
     */
    public static <T extends Comparable<T>> T max(T x, T y) {      return (x.compareTo(y) > 0) ? x : y;    }
    /**
     * return min value of any comparable type
     */
    public static <T extends Comparable<T>> T min(T x, T y) {      return (x.compareTo(y) < 0) ? x : y;    }
   
    /**
     * return max value of any comparable type of 3 values
     */
    public static <T extends Comparable<T>> T max(T x, T y, T z) {    	return max(max(x,y),z);    }
    /**
     * return min value of any comparable type
     */
    public static <T extends Comparable<T>> T min(T x, T y, T z) {    	return min(min(x,y),z);     }
    
	/**
	 * Range checking of value, as double
	 * @param x x value to check
	 * @param min minimum x value (inclusive)
	 * @param max maximum x value (inclusive)
	 * @return if value is in specified range
	 */
    public static boolean inRange(double x, double min, double max) {return (x >= min)&&(x <= max); }	
    
	/**
	 * Range checking of value, as floats
	 * @param x x value to check
	 * @param min minimum x value (inclusive)
	 * @param max maximum x value (inclusive)
	 * @return if value is in specified range
	 */
    public static boolean inRange(float x, float min, float max) {return (x >= min)&&(x <= max); }	
    
	/**
	 * Range checking of value, as integers
	 * @param x x value to check
	 * @param min minimum x value (inclusive)
	 * @param max maximum x value (inclusive)
	 * @return if value is in specified range
	 */
    public static boolean inRange(int x, int min, int max) {return (x >= min)&&(x <= max); }	
    
	/**
	 * Range checking of value, as longs
	 * @param x x value to check
	 * @param min minimum x value (inclusive)
	 * @param max maximum x value (inclusive)
	 * @return if value is in specified range
	 */
    public static boolean inRange(long x, long min, long max) {return (x >= min)&&(x <= max); }	
    
	/**
	 * Range checking of value of any types that extend Comparable interface
	 * @param x x value to check
	 * @param min minimum x value (inclusive)
	 * @param max maximum x value (inclusive)
	 * @return if value is in specified range
	 */
    public static <T extends Comparable<T>> boolean inRange(T x, T min, T max) {return (x.compareTo(min) >=0) && (x.compareTo(max) <= 0); }	
       
	/**
	 * 2d range checking of point, represented as doubles
	 * @param x x value to check
	 * @param y y value to check
	 * @param minX minimum x value (inclusive)
	 * @param minY minimum y value (inclusive)
	 * @param maxX maximum x value (inclusive)
	 * @param maxY maximum y value (inclusive)
	 * @return if both values are in specified range
	 */
    public static boolean ptInRange(double x, double y, double minX, double minY, double maxX, double maxY) {return ((x >= minX)&&(x <= maxX)&&(y >= minY)&&(y <= maxY)); }	
    
	/**
	 * 2d range checking of point, represented as floats
	 * @param x x value to check
	 * @param y y value to check
	 * @param minX minimum x value (inclusive)
	 * @param minY minimum y value (inclusive)
	 * @param maxX maximum x value (inclusive)
	 * @param maxY maximum y value (inclusive)
	 * @return if both values are in specified range
	 */
    public static boolean ptInRange(float x, float y, float minX, float minY, float maxX, float maxY) {return ((x >= minX)&&(x <= maxX)&&(y >= minY)&&(y <= maxY)); }	
    
	/**
	 * 2d range checking of point, represented as integers
	 * @param x x value to check
	 * @param y y value to check
	 * @param minX minimum x value (inclusive)
	 * @param minY minimum y value (inclusive)
	 * @param maxX maximum x value (inclusive)
	 * @param maxY maximum y value (inclusive)
	 * @return if both values are in specified range
	 */
    public static boolean ptInRange(int x, int y, int minX, int minY, int maxX, int maxY) {return ((x >= minX)&&(x <= maxX)&&(y >= minY)&&(y <= maxY)); }	
    
	/**
	 * 2d range checking of point, represented as longs
	 * @param x x value to check
	 * @param y y value to check
	 * @param minX minimum x value (inclusive)
	 * @param minY minimum y value (inclusive)
	 * @param maxX maximum x value (inclusive)
	 * @param maxY maximum y value (inclusive)
	 * @return if both values are in specified range
	 */
    public static boolean ptInRange(long x, long y, long minX, long minY, long maxX, long maxY) {return ((x >= minX)&&(x <= maxX)&&(y >= minY)&&(y <= maxY)); }	
    
	/**
	 * 2d range checking of points of any types that extend Comparable interface
	 * @param x x value to check
	 * @param y y value to check
	 * @param minX minimum x value (inclusive)
	 * @param minY minimum y value (inclusive)
	 * @param maxX maximum x value (inclusive)
	 * @param maxY maximum y value (inclusive)
	 * @return if both values are in specified range
	 */
    public static <T extends Comparable<T>> boolean ptInRange(T x, T y, T minX, T minY, T maxX, T maxY) {
    	return ((x.compareTo(minX) >=0) && (x.compareTo(maxX) <= 0) && (y.compareTo(minY) >=0) && (y.compareTo(maxY) <= 0) );
    }	
	
	/**
	 * n-dim range checking of array of values of doubles.
	 * @param vals array of vals to check
	 * @param min array of min bounds (inclusive)
	 * @param max array of max bounds (inclusive)
	 * @return whether entire value array is within bounds
	 */
	public static boolean valuesInRange(double[] vals, double[] mins, double[] maxs) {
		if((vals.length != mins.length) || (vals.length != maxs.length)){return false;}	//insufficient bounds passed
		int len = vals.length;
		for(int i=0;i<len;++i) {
			if((vals[i] <mins[i]) || (vals[i] > maxs[i])){	return false;}		
		}		
		return true;
	}
	
	/**
	 * n-dim range checking of array of values of floats.
	 * @param vals array of vals to check
	 * @param min array of min bounds (inclusive)
	 * @param max array of max bounds (inclusive)
	 * @return whether entire value array is within bounds
	 */
	public static boolean valuesInRange(float[] vals, float[] mins, float[] maxs) {
		if((vals.length != mins.length) || (vals.length != maxs.length)){return false;}	//insufficient bounds passed
		int len = vals.length;
		for(int i=0;i<len;++i) {
			if((vals[i] <mins[i]) || (vals[i] > maxs[i])){	return false;}		
		}		
		return true;
	}
	
	/**
	 * n-dim range checking of array of values of integers.
	 * @param vals array of vals to check
	 * @param min array of min bounds (inclusive)
	 * @param max array of max bounds (inclusive)
	 * @return whether entire value array is within bounds
	 */
	public static boolean valuesInRange(int[] vals, int[] mins, int[] maxs) {
		if((vals.length != mins.length) || (vals.length != maxs.length)){return false;}	//insufficient bounds passed
		int len = vals.length;
		for(int i=0;i<len;++i) {
			if((vals[i] <mins[i]) || (vals[i] > maxs[i])){	return false;}		
		}		
		return true;
	}
	
	/**
	 * n-dim range checking of array of values of longs.
	 * @param vals array of vals to check
	 * @param min array of min bounds (inclusive)
	 * @param max array of max bounds (inclusive)
	 * @return whether entire value array is within bounds
	 */
	public static boolean valuesInRange(long[] vals, long[] mins, long[] maxs) {
		if((vals.length != mins.length) || (vals.length != maxs.length)){return false;}	//insufficient bounds passed
		int len = vals.length;
		for(int i=0;i<len;++i) {
			if((vals[i] < mins[i]) || (vals[i] > maxs[i])){	return false;}		
		}		
		return true;
	}
	
	/**
	 * n-dim range checking of values of types that extend Comparable interface
	 * @param vals array of vals to check
	 * @param min array of min bounds (inclusive)
	 * @param max array of max bounds (inclusive)
	 * @return whether entire value array is within bounds
	 */
	public static <T extends Comparable<T>> boolean valuesInRange(T[] vals, T[] mins, T[] maxs) {
		if((vals.length != mins.length) || (vals.length != maxs.length)){return false;}	//insufficient bounds passed
		int len = vals.length;
		for(int i=0;i<len;++i) {
			if((vals[i].compareTo(mins[i]) < 0)  || (vals[i].compareTo(maxs[i]) > 0)){	return false;}		
		}		
		return true;
	}
     
}//math utils

