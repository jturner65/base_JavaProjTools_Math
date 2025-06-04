package base_Math_Objects;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

import base_Math_Objects.vectorObjs.doubles.myPoint;
import base_Math_Objects.vectorObjs.doubles.myVector;
import base_Math_Objects.vectorObjs.floats.myPointf;
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
	public static double
		THIRD = 1.0/3.0,
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
	public static float 
		THIRD_F = 1.0f/3.0f,
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
    protected static int MAX_DIGITS_EXP = 677;
    protected static int MAX_DIGITS_10 = 294; // ~ MAX_DIGITS_EXP/LN(10)
    protected static int MAX_DIGITS_2 = 977; // ~ MAX_DIGITS_EXP/LN(2)
        
    /**
     * # of pre-computed cosine and sine values.  
     */
    public static float numTrigVals = 36.0f;
    /**
     * Pre-computed array of {@value #numTrigVals} cosine values 0->2pi
     */
    public static double[] preCalcCosVals;
    /**
     * Pre-computed array of {@value #numTrigVals} sine values 0->2pi
     */
    public static double[] preCalcSinVals;
    
    private static float deltaThet = TWO_PI_F/numTrigVals,
			finalThet = TWO_PI_F+deltaThet;
    /**
     * Pre-computed array of {@value #numTrigVals} float cosine values 0->2pi
     */
    public static float[] preCalcCosVals_f;
    /**
     * Pre-computed array of {@value #numTrigVals} float sine values 0->2pi
     */
    public static float[] preCalcSinVals_f;
    
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
	
	//public static int O_FWD = 0, O_RHT = 1,  O_UP = 2;
	/**
	 * build axis angle orientation from passed orientation matrix
	 * @param orientation array of 3 vectors corresponding to orientation vectors
	 * @param O_FWD idx of forward orientation
	 * @param O_RHT idx of right orientation
	 * @param O_UP idx of up orientation
	 * @return axis-angle representation of orientation
	 */
	public static float[] toAxisAngle(myVectorf[] orientation, int O_FWD, int O_RHT, int O_UP) {
		float angle,s,x=INV_SQRT_2_F,y=INV_SQRT_2_F,z=INV_SQRT_2_F;
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
	public static myVector normToPlane(myPoint A, myPoint B, myPoint C) {
		return myVector._cross(new myVector(A,B),new myVector(A,C)); 
	};   // normal to triangle (A,B,C), not normalized (proportional to area)
	
	/**
	 * Calculate normal to planed described by triangle ABC, non-normalized (proportional to area)
	 * @param A, B, C verts of triangle
	 * @return
	 */
	public static myVectorf normToPlane(myPointf A, myPointf B, myPointf C) {
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
	public static int floor(float x){	return x>=0 ? (int)x : (int)x-1;} 
	 
	/**
	 * Floor of passed double
	 * @param x
	 * @return
	 */
	public static int floor(double x){ return x>=0 ? (int)x : (int)x-1;} 
	
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
     
	
	/**
	 * point normal form of plane - give a normal and a point and receive a 4-element array of the coefficients of the plane equation.
	 * @param _n unit normal of plane
	 * @param _p point on plane
	 * @return eq : coefficients of the plane equation in the form eq[0]*x + eq[1]*y + eq[2]*z + eq[3] = 0
	 */
	public static float[] getPlanarEqFromPointAndNorm(myVector _n, myPoint _p) {
		return new float[] { (float)_n.x,(float) _n.y, (float)_n.z, (float) (_n.x * -_p.x + _n.y * -_p.y + _n.z * -_p.z)};
	}
	

	/**
	 * point normal form of plane - give a normal and a point and receive a 4-element array of the coefficients of the plane equation.
	 * @param _n unit normal of plane
	 * @param _p point on plane
	 * @return eq : coefficients of the plane equation in the form eq[0]*x + eq[1]*y + eq[2]*z + eq[3] = 0
	 */
	public static float[] getPlanarEqFromPointAndNorm(myVectorf _n, myPointf _p) {
		return new float[] { _n.x, _n.y, _n.z, (_n.x * -_p.x + _n.y * -_p.y + _n.z * -_p.z)};
	}	
	

	/**
	 * build a frame based on passed normal given two passed points
	 * @param A
	 * @param B
	 * @param I normal to build frame around
	 * @return vec array of {AB, Normal, Tangent}
	 */
	public static myVector[] buildFrameAroundNormal(myPoint A, myPoint B, myVector C) {
		myVector V = new myVector(A,B);
		myVector norm = myVector._findNormToPlane(C,V);		
		myVector tan = norm._cross(V)._normalize(); 
		return new myVector[] {V,norm,tan};		
	}
	
	/**
	 * build a frame based on passed normal given two passed points
	 * @param A
	 * @param B
	 * @param I normal to build frame around
	 * @return vec array of {AB, Normal, Tangent}
	 */
	public static myVectorf[] buildFrameAroundNormal(myPointf A, myPointf B, myVectorf C) {
		myVectorf V = new myVectorf(A,B);
		myVectorf norm = myVectorf._findNormToPlane(C,V);
		myVectorf tan = norm._cross(V)._normalize(); 
		return new myVectorf[] {V,norm,tan};		
	}
	
	/**
	 * Derive the points of a cylinder of radius r around axis through A and B
	 * @param A center point of endcap
	 * @param B center point of endcap
	 * @param r desired radius of cylinder
	 * @return array of points for cylinder
	 */
	public static myPoint[] buildCylVerts(myPoint A, myPoint B, double r, myVector norm) {
		myVector[] frame = buildFrameAroundNormal(A, B, norm);
		myPoint[] resList = new myPoint[2 * preCalcCosVals.length];
		double rca, rsa;
		int idx = 0;
		for(int i = 0; i<preCalcCosVals.length; ++i) {
			rca = r*preCalcCosVals[i];
			rsa = r*preCalcCosVals[i];
			resList[idx++] = myPoint._add(A,rca,frame[1],rsa,frame[2]); 
			resList[idx++] = myPoint._add(A,rca,frame[1],rsa,frame[2],1,frame[0]);				
		}
		return resList;
	}//build list of all cylinder vertices 
	
	/**
	 * Derive the points of a cylinder of radius r around axis through A and B
	 * @param A center point of endcap
	 * @param B center point of endcap
	 * @param r desired radius of cylinder
	 * @return array of points for cylinder
	 */
	public static myPointf[] buildCylVerts(myPointf A, myPointf B, float r, myVectorf norm) {
		myVectorf[] frame = buildFrameAroundNormal(A, B, norm);
		myPointf[] resList = new myPointf[2 * preCalcCosVals.length];
		float rca, rsa;
		int idx = 0;
		for(int i = 0; i<preCalcCosVals.length; ++i) {
			rca = r*preCalcCosVals_f[i];
			rsa = r*preCalcSinVals_f[i];
			resList[idx++] = myPointf._add(A,rca,frame[1],rsa,frame[2]); 
			resList[idx++] = myPointf._add(A,rca,frame[1],rsa,frame[2],1,frame[0]);				
		}	
		return resList;
	}//build list of all cylinder vertices 
	
	/**
	 * Build a set of n points inscribed on a circle centered at p in plane I,J
	 * @param p center point
	 * @param r circle radius
	 * @param I, J axes of plane
	 * @param n # of points
	 * @return array of n equal-arc-length points centered around p
	 */
	public static myPoint[] buildCircleInscribedPoints(myPoint p, double r, myVector I, myVector J, int n) {
		myPoint[] pts = new myPoint[n];
		pts[0] = new myPoint(p,r,myVector._unit(I));
		double a = (TWO_PI)/(1.0*n); 
		for(int i=1;i<n;++i){pts[i] = pts[i-1].rotMeAroundPt(a,J,I,p);}
		return pts;
	}
	/**
	 * Build a set of n points inscribed on a circle centered at p in plane I,J
	 * @param p center point
	 * @param r circle radius
	 * @param I, J axes of plane
	 * @param n # of points
	 * @return array of n equal-arc-length points centered around p
	 */
	public static myPointf[] buildCircleInscribedPoints(myPointf p, float r, myVectorf I, myVectorf J, int n) {
		myPointf[] pts = new myPointf[n];
		pts[0] = new myPointf(p,r,myVectorf._unit(I));
		float a = (TWO_PI_F)/(1.0f*n);
		for(int i=1;i<n;++i){pts[i] = pts[i-1].rotMeAroundPt(a,J,I,p);}
		return pts;
	}
	
	
	/**
	 * 
	 * @param A
	 * @param B
	 * @param C
	 * @param t
	 * @return
	 */
	public static myPoint PtOnSpiral(myPoint A, myPoint B, myPoint C, double t) {
		//center is coplanar to A and B, and coplanar to B and C, but not necessarily coplanar to A, B and C
		//so center will be coplanar to mp(A,B) and mp(B,C) - use mpCA midpoint to determine plane mpAB-mpBC plane?
		myPoint mAB = new myPoint(A,.5, B), mBC = new myPoint(B,.5, C), mCA = new myPoint(C,.5, A);
		myVector mI = myVector._unit(mCA,mAB), mTmp = myVector._cross(mI,myVector._unit(mCA,mBC)), mJ = myVector._unit(mTmp._cross(mI));	//I and J are orthonormal
		double a =spiralAngle(A,B,B,C), s =spiralScale(A,B,B,C);
		
		//myPoint G = spiralCenter(a, s, A, B, mI, mJ); 
		myPoint G = spiralCenter(A, mAB, B, mBC); 
		//return new myPoint(G, Math.pow(s,t), R(A,t*a,mI,mJ,G));
		return new myPoint(G, Math.pow(s,t), A.rotMeAroundPt(t*a,mI,mJ,G));
	}	
	/**
	 * 
	 * @param A
	 * @param B
	 * @param C
	 * @param D
	 * @return
	 */
	public static double spiralAngle(myPoint A, myPoint B, myPoint C, myPoint D) {return myVector._angleBetween(new myVector(A,B),new myVector(C,D));}
	/**
	 * 
	 * @param A
	 * @param B
	 * @param C
	 * @param D
	 * @return
	 */
	public static double spiralScale(myPoint A, myPoint B, myPoint C, myPoint D) {return myPoint._dist(C,D)/ myPoint._dist(A,B);}

	/**
	 * spiral given 4 points, AB and CD are edges corresponding through rotation
	 * @param A
	 * @param B
	 * @param C
	 * @param D
	 * @return
	 */
	public static myPoint spiralCenter(myPoint A, myPoint B, myPoint C, myPoint D) {         // new spiral center
		myVector AB=new myVector(A,B), CD=new myVector(C,D), AC=new myVector(A,C);
		double m=CD.magn/AB.magn, n=CD.magn*AB.magn;		
		myVector rotAxis = myVector._unit(AB._cross(CD));		//expect ab and ac to be coplanar - this is the axis to rotate around to find flock
		
		myVector rAB = myVector._rotAroundAxis(AB, rotAxis, HALF_PI_F);
		double c=AB._dot(CD)/n,	s=rAB._dot(CD)/n;
		double AB2 = AB._dot(AB), a=AB._dot(AC)/AB2, b=rAB._dot(AC)/AB2, x=(a-m*( a*c+b*s)), y=(b-m*(-a*s+b*c)), d=1+m*(m-2*c);  if((c!=1)&&(m!=1)) { x/=d; y/=d; };
		return new myPoint(new myPoint(A,x,AB),y,rAB);
	}
	/**
	 * 
	 * @param A
	 * @param B
	 * @param C
	 * @param t
	 * @return
	 */
	public static myPointf PtOnSpiral(myPointf A, myPointf B, myPointf C, float t) {
		//center is coplanar to A and B, and coplanar to B and C, but not necessarily coplanar to A, B and C
		//so center will be coplanar to mp(A,B) and mp(B,C) - use mpCA midpoint to determine plane mpAB-mpBC plane?
		myPointf mAB = new myPointf(A,.5f, B), mBC = new myPointf(B,.5f, C), mCA = new myPointf(C,.5f, A);
		myVectorf mI = myVectorf._unit(mCA,mAB), mTmp = myVectorf._cross(mI,myVectorf._unit(mCA,mBC)), mJ = myVectorf._unit(mTmp._cross(mI));	//I and J are orthonormal
		float a =spiralAngle(A,B,B,C), s =spiralScale(A,B,B,C);
		
		//myPoint G = spiralCenter(a, s, A, B, mI, mJ); 
		myPointf G = spiralCenter(A, mAB, B, mBC); 
		return new myPointf(G, (float)Math.pow(s,t), A.rotMeAroundPt(t*a,mI,mJ,G));
	}	
	/**
	 * 
	 * @param A
	 * @param B
	 * @param C
	 * @param D
	 * @return
	 */
	public static float spiralAngle(myPointf A, myPointf B, myPointf C, myPointf D) {return myVectorf._angleBetween(new myVectorf(A,B),new myVectorf(C,D));}
	/**
	 * 
	 * @param A
	 * @param B
	 * @param C
	 * @param D
	 * @return
	 */
	public static float spiralScale(myPointf A, myPointf B, myPointf C, myPointf D) {return myPointf._dist(C,D)/ myPointf._dist(A,B);}
	
	/**
	 * spiral given 4 points, AB and CD are edges corresponding through rotation
	 * @param A
	 * @param B
	 * @param C
	 * @param D
	 * @return
	 */
	public static myPointf spiralCenter(myPointf A, myPointf B, myPointf C, myPointf D) {         // new spiral center
		myVectorf AB=new myVectorf(A,B), CD=new myVectorf(C,D), AC=new myVectorf(A,C);
		float m=CD.magn/AB.magn, n=CD.magn*AB.magn;		
		myVectorf rotAxis = myVectorf._unit(AB._cross(CD));		//expect ab and ac to be coplanar - this is the axis to rotate around to find flock
		
		myVectorf rAB = myVectorf._rotAroundAxis(AB, rotAxis, HALF_PI_F);
		float c=AB._dot(CD)/n,	s=rAB._dot(CD)/n;
		float AB2 = AB._dot(AB), a=AB._dot(AC)/AB2, b=rAB._dot(AC)/AB2, x=(a-m*( a*c+b*s)), y=(b-m*(-a*s+b*c)), d=1+m*(m-2*c);  if((c!=1)&&(m!=1)) { x/=d; y/=d; };
		return new myPointf(new myPointf(A,x,AB),y,rAB);
	}

	/**
	 * Return intersection point in plane described by ABC of vector T through point E
	 * @param E point within cast vector/ray
	 * @param T directional vector/ray
	 * @param A point describing plane
	 * @param B point describing plane
	 * @param C point describing plane
	 * @return
	 */
	public static myPoint intersectPlane(myPoint E, myVectorf T, myPointf A, myPointf B, myPointf C) {
		myPointf res = intersectPlane(new myPointf(E.x, E.y, E.z), T, A, B, C);
		return new myPoint(res.x, res.y, res.z);
	}//intersectPlane
	/**
	 * Return intersection point in plane described by ABC of vector T through point E
	 * @param E point within cast vector/ray
	 * @param T directional vector/ray
	 * @param A point describing plane
	 * @param B point describing plane
	 * @param C point describing plane
	 * @return
	 */
	public static myPointf intersectPlane(myPointf E, myVectorf T, myPointf A, myPointf B, myPointf C) {
		//vector through point and planar point
		myVectorf EA=new myVectorf(E,A); 
		//planar vectors
		myVectorf AB=new myVectorf(A,B), AC=new myVectorf(A,C);
		//find planar norm
		myVectorf ACB = AC._cross(AB);
		//project 
		float t = (EA._dot(ACB) / T._dot(ACB));		
		return (myPointf._add(E,t,T));		
	}//intersectPlane
	/**
	 * Return intersection point in plane described by ABC of vector T through point E
	 * @param E point within cast vector/ray
	 * @param T directional vector/ray
	 * @param A point describing plane
	 * @param B point describing plane
	 * @param C point describing plane
	 * @return
	 */
	public static myPoint intersectPlane(myPoint E, myVector T, myPoint A, myPoint B, myPoint C) {
		//vector through point and planar point
		myVector EA=new myVector(E,A); 
		//planar vectors
		myVector AB=new myVector(A,B), AC=new myVector(A,C);
		//find planar norm
		myVector ACB = AC._cross(AB);
		//project 
		double t = (EA._dot(ACB) / T._dot(ACB));		
		return (myPoint._add(E,t,T));		
	}//intersectPlane		
	/**
	 * if ray from E along V intersects sphere at C with radius r, return t when intersection occurs
	 * @param E
	 * @param V
	 * @param C
	 * @param r
	 * @return t value along vector V where first intersection occurs
	 */
	public static double intersectPt(myPoint E, myVector V, myPoint C, double r) { 
		myVector Vce = new myVector(C,E);
		double ta = 2 * V._dot(V),
				b = 2 * V._dot(Vce), 
				c = Vce._dot(Vce) - (r*r),
				radical = (b*b) - 2 *(ta) * c;		//b^2 - 4ac
		if(radical < 0) return -1;
		double sqrtRad = Math.sqrt(radical);
		double t1 = (b + sqrtRad)/ta, t2 = (b - sqrtRad)/ta;
		if (t1 < t2) {return t1 > 0 ? t1 : t2;}	
		return t2 > 0 ? t2 : t1;	
		//return ((t1 > 0) && (t2 > 0) ? min(t1, t2) : ((t1 < 0 ) ? ((t2 < 0 ) ?-1 : t2) : t1) );
	}	
	
	////////////////////////////////////////////
	/// Synchronized Random functions
	/**
	 * Get random next double between 0.0 (inclusive) and 1.0 (exclusive). Synchronized.
	 * @return
	 */
	public static synchronized double randomDouble() {return ThreadLocalRandom.current().nextDouble();}	

	/**
	 * Get random next double between 0 and @param high. Synchronized.
	 * @param high upper (exclusionary) bound for random draw
	 * @return
	 */
	public static synchronized double randomDouble(double high) {return ThreadLocalRandom.current().nextDouble(high);}	

	/**
	 * Get random next double between @param low and @param high. Synchronized.
	 * @param low lower, inclusive, bound for random draw
	 * @param high upper (exclusionary) bound for random draw
	 * @return
	 */
	public static synchronized double randomDouble(double low, double high) {return ThreadLocalRandom.current().nextDouble(low, high);}

	/**
	 * Get random next float between 0.0f (inclusive) and 1.0f (exclusive). Synchronized.
	 * @return
	 */
	public static synchronized float randomFloat() {return ThreadLocalRandom.current().nextFloat();}

	/**
	 * Get random next float between 0 and @param high. Synchronized.
	 * @param high upper (exclusionary) bound for random draw
	 * @return
	 */
	public static synchronized float randomFloat(float high) {return ThreadLocalRandom.current().nextFloat(high);}

	/**
	 * Get random next float between @param low and @param high. Synchronized.
	 * @param low lower, inclusive, bound for random draw
	 * @param high upper (exclusionary) bound for random draw
	 * @return
	 */
	public static synchronized float randomFloat(float low, float high) {return ThreadLocalRandom.current().nextFloat(low, high);}
	
	/**
	 * Get random next int between -Integer.MAX_VALUE and +Integer.MAX_VALUE. Synchronized.
	 * @return
	 */
	public static synchronized int randomInt() {return ThreadLocalRandom.current().nextInt();}
	
	/**
	 * Get random next int between 0 and @param high. Synchronized.
	 * @param high upper (exclusionary) bound for random draw
	 * @return
	 */
	public static synchronized int randomInt(int high) {return ThreadLocalRandom.current().nextInt(high);}
	
	/**
	 * Get random next int between @param low and @param high. Synchronized.
	 * @param low lower, inclusive, bound for random draw
	 * @param high upper (exclusionary) bound for random draw
	 * @return
	 */
	public static synchronized int randomInt(int low, int high) {return ThreadLocalRandom.current().nextInt(low, high);}
	
	/**
	 * Get random next long between -Long.MAX_VALUE and +Long.MAX_VALUE. Synchronized. 
	 * Due to only using a 48-bit seed, not all possible long values are generated, but those 
	 * that are generated are stated to be uniformly distributed in the given range.
	 * @return
	 */
	public static synchronized long randomLong() {return ThreadLocalRandom.current().nextLong();}
	
	/**
	 * Get random next long between 0 and @param high. Synchronized.
	 * @param high upper (exclusionary) bound for random draw
	 * @return
	 */
	public static synchronized long randomLong(long high) {return ThreadLocalRandom.current().nextLong(high);}
	
	/**
	 * Get random next long between @param low and @param high. Synchronized.
	 * @param low lower, inclusive, bound for random draw
	 * @param high upper (exclusionary) bound for random draw
	 * @return
	 */
	public static synchronized long randomLong(long low, long high) {return ThreadLocalRandom.current().nextLong(low, high);}
	
	/**
	 * Get a random uniform double between 0 and 1.
	 * @return
	 */
	public static synchronized double randomUniform01() {
		//-Integer.MAX_VALUE to Integer.MAX_VALUE
    	int val = ThreadLocalRandom.current().nextInt();
    	return .5+ .5 * val/Integer.MAX_VALUE;
	}
	
	/**
	 * Random int between @low and @high that is a legal color channel (i.e. [0,255])
	 * @param low lower, inclusionary bound of red channel
	 * @param high
	 * @return
	 */
	public static final int randomIntClrChan(int low, int high) {
		low = (low < 0 ? 0 : low);
		high = (high > 256 ? 256 : high);
		return randomInt(low,high);
	}
	/**
	 * Return a random 4-element int array of values suitable for color with each color between minX, inclusive, and maxX, 
	 * exclusive (ie. use 256 for maxX to get highest bound), with given alpha value.
	 * @param minR lower, inclusionary bound of red channel
	 * @param maxR upper, exclusionary, bound of red channel
	 * @param minG lower, inclusionary bound of green channel
	 * @param maxG upper, exclusionary, bound of green channel
	 * @param minB lower, inclusionary bound of blue channel
	 * @param maxB upper, exclusionary, bound of blue channel
	 * @param alpha alpha channel value to use
	 * @return
	 */
	public static final int[] randomIntClrAra(int minR, int maxR, int minG, int maxG, int minB, int maxB, int alpha) {
		return new int[] {randomIntClrChan(minR, maxR), randomIntClrChan(minG,maxG), randomIntClrChan(minB,maxB), alpha};
	}	

	/**
	 * Return a random 4-element int array of values suitable for color with each channel less than 
	 * the passed upper bound (ie. use 256 for maxX to get highest bound), with given alpha value.
	 * @param maxR upper, exclusionary, bound of red channel
	 * @param maxG upper, exclusionary, bound of green channel
	 * @param maxB upper, exclusionary, bound of blue channel
	 * @param alpha alpha channel value to use
	 * @return
	 */
	public static final int[] randomIntClrAra(int maxR, int maxG, int maxB, int alpha) {return randomIntClrAra(0, maxR, 0, maxG, 0, maxB, alpha);}
	
	/**
	 * Return a random 4-element int array of values suitable for color where each channel is between [0,255], with given alpha value
	 * @param alpha alpha channel value to use
	 * @return
	 */
	public static final int[] randomIntClrAra(int alpha) {return randomIntClrAra(0, 256, 0, 256, 0, 256, alpha);}	
		
	/**
	 * Return a random 4-element int array of values suitable for color with each color between minX, inclusive, and maxX, 
	 * exclusive (ie. use 256 for maxX to get highest bound) with alpha == 255
	 * @param minR lower, inclusionary bound of red channel
	 * @param maxR upper, exclusionary, bound of red channel
	 * @param minG lower, inclusionary bound of green channel
	 * @param maxG upper, exclusionary, bound of green channel
	 * @param minB lower, inclusionary bound of blue channel
	 * @param maxB upper, exclusionary, bound of blue channel
	 * @return
	 */
	public static final int[] randomIntClrAra(int minR, int maxR, int minG, int maxG, int minB, int maxB) {return randomIntClrAra(minR, maxR,minG,maxG,minB,maxB, 255);}

	/**
	 * Return a random 4-element int array of values suitable for color with each channel less than 
	 * the passed upper bound (ie. use 256 for maxX to get highest bound) with alpha == 255
	 * @param maxR upper, exclusionary, bound of red channel
	 * @param maxG upper, exclusionary, bound of green channel
	 * @param maxB upper, exclusionary, bound of blue channel
	 * @return
	 */
	public static final int[] randomIntClrAra(int maxR, int maxG, int maxB) {return randomIntClrAra(0, maxR, 0, maxG, 0, maxB, 255);}

	/**
	 * Return a random 4-element int array of values suitable for color where each channel is between [0,255], with alpha == 255
	 * @return
	 */
	public static final int[] randomIntClrAra() {return randomIntClrAra(0, 256, 0, 256, 0, 256, 255);}
	
	/**
	 * Return a random 4-element int array of values suitable for color with each channel between [(r:75,g:30,b:80),255], with alpha == 255
	 * @return
	 */
	public static final int[] randomIntBrightClrAra() {	return randomIntClrAra(75, 256, 30, 256, 80, 256, 255);}

	////////////////////////////////////////////
	/// End Synchronized Random functions	
	
	/**
	 * Find a random position within a sphere centered at 0 of radius rad, using spherical coords as rand axes
	 * @param rad
	 * @return
	 */
	public static myPointf getRandPosInSphere(double rad){ return getRandPosInSphere(rad, new myPointf());}
	/**
	 * Find a random position within a sphere centered at ctr of radius rad, using spherical coords as rand axes
	 * @param rad
	 * @param ctr
	 * @return
	 */
	public static myPointf getRandPosInSphere(double rad, myPointf ctr){
		myPointf pos = new myPointf();
		double u = randomDouble(0,1),	
			cosTheta = randomDouble(-1,1),
			phi = randomDouble(0,TWO_PI_F),
			r = rad * Math.pow(u, THIRD_F),
			rSinTheta = r * (Math.sqrt(1.0 - (cosTheta * cosTheta)));			
		pos.set(rSinTheta * Math.cos(phi), rSinTheta * Math.sin(phi),cosTheta*r);
		pos._add(ctr);
		return pos;
	}
	/**
	 * Find a random position on a sphere's surface centered at 0 of radius rad, using spherical coords as rand axes
	 * @param rad
	 * @return
	 */
	public static myPointf getRandPosOnSphere(double rad){ return getRandPosOnSphere(rad, new myPointf());}
	/**
	 * Find a random position on a sphere's surface centered at ctr of radius rad, using spherical coords as rand axes
	 * @param rad
	 * @param ctr
	 * @return
	 */
	public static myPointf getRandPosOnSphere(double rad, myPointf ctr){
		myPointf pos = new myPointf();
		double 	cosTheta = randomDouble(-1,1),
				phi = randomDouble(0,TWO_PI_F), 
				rSinTheta = rad* (Math.sqrt(1.0 - (cosTheta * cosTheta)));
		pos.set(rSinTheta * Math.cos(phi), rSinTheta * Math.sin(phi),cosTheta * rad);
		pos._add(ctr);
		return pos;
	}
	/**
	 * Find a random position on a sphere's surface centered at 0 of radius rad, using spherical coords as rand axes
	 * @param rad
	 * @return
	 */
	public static myPoint getRandPosOnSphereDouble(double rad){ return getRandPosOnSphereDouble(rad, new myPoint());}
	/**
	 * Find a random position on a sphere's surface centered at ctr of radius rad, using spherical coords as rand axes
	 * @param rad
	 * @param ctr
	 * @return
	 */
	public static myPoint getRandPosOnSphereDouble(double rad, myPoint ctr){
		myPoint pos = new myPoint();
		double 	cosTheta = randomDouble(-1,1),
				phi = randomDouble(0,TWO_PI_F), 
				rSinTheta = rad* (Math.sqrt(1.0 - (cosTheta * cosTheta)));
		pos.set(rSinTheta * Math.cos(phi), rSinTheta * Math.sin(phi),cosTheta * rad);
		pos._add(ctr);
		return pos;
	}
	/** 
	 * convert from spherical coords to cartesian. Returns Array :
	 * 	idx0 : norm of vector through point from origin
	 * 	idx1 : point
	 * @param rad
	 * @param thet
	 * @param phi
	 * @param scaleZ scaling factor to make ellipsoid
	 * @return ara : norm, surface point == x,y,z of coords passed
	 */
	public static myVectorf[] getXYZFromRThetPhi(double rad, double thet, double phi, double scaleZ) {
		double sinThet = Math.sin(thet);	
		myVectorf[] res = new myVectorf[2];
		res[1] = new myVectorf(sinThet * Math.cos(phi) * rad, sinThet * Math.sin(phi) * rad,Math.cos(thet)*rad*scaleZ);
		res[0] = myVectorf._normalize(res[1]);
		return res;
	}//	
	
	/** 
	 * builds a list of N regularly placed vertices and normals for a sphere of radius rad centered at ctr
	 * @param rad radius of sphere
	 * @param N # of verts we want in result
	 * @param scaleZ scaling factor for ellipsoid
	 * @return list of points (as vectors) where each entry is a tuple of norm/point
	 */
	public static myVectorf[][] getRegularSphereList(float rad, int N, float scaleZ) {
		ArrayList<myVectorf[]> res = new ArrayList<myVectorf[]>();
		//choose 1 point per dArea, where dArea is area of sphere parsed into N equal portions
		double lclA = 4*PI/N, lclD = Math.sqrt(lclA);
		int Mthet = (int) Math.round(PI/lclD), Mphi;
		double dThet = PI/Mthet, dPhi = lclA/dThet, thet, phi, twoPiOvDPhi = TWO_PI/dPhi;
		for(int i=0;i<Mthet;++i) {
			thet = dThet * (i + 0.5f);
			Mphi = (int) Math.round(twoPiOvDPhi * Math.sin(thet));
			for (int j=0;j<Mphi; ++j) { 
				phi = (TWO_PI*j)/Mphi;		
				res.add(getXYZFromRThetPhi(rad, thet, phi, scaleZ));
			}
		}
		return res.toArray(new myVectorf[0][]);
	}//getRegularSphereList	

	
}//math utils

