package base_Math_Objects.vectorObjs.tuples;


//only works for comparables
public class Tuple<X,Y> implements Comparable<Tuple<X,Y>> { 
	//private static final int[] seeds = new int[] {7919, 1699};
    public final X x;    public final Y y;  private final Float sqmag; private final int hash;
    public Tuple(X x, Y y) {    this.x = x;   this.y = y; sqmag = calcSqMag(); hash= this.hashCode(); }
    public Tuple(Tuple<X,Y> _t) {   this( _t.x,_t.y);  }
    public String toCSVString() {	return "("+x+"|"+y+")";}
    public String toString() {      return "(" + x + "," + y + ")";  }
    @SuppressWarnings("unchecked")
	public boolean equals(Object _o) {  if (_o == null) {return false;} if (_o == this) { return true; } if (!(_o instanceof Tuple)){ return false; } Tuple<X,Y> o = (Tuple<X,Y>) _o;  return o.x.equals(this.x) && o.y.equals(this.y);  }
    public int hashCode() { 
//    	Random random = new Random(seeds[0] + x.hashCode());
//    	long result = random.nextInt();
//    	random.setSeed(seeds[1] + y.hashCode());
//    	result += random.nextInt();
//    	return (int) (result % Integer.MAX_VALUE);
    	
    	int result = 97 + ((x == null) ? 0 : x.hashCode());return 97 * result + ((y == null) ? 0 : y.hashCode());      	
    }
    public Float getSqMag() {return sqmag;}
	public Float calcSqMag(){if((x != null) && (y != null)) { return 1.0f*((x.hashCode()*x.hashCode()) + (y.hashCode()*y.hashCode()));} else {return null;}}
 	@Override
	public int compareTo(Tuple<X, Y> arg0) {//not a good measure - need to first use dist 		
 		if (this.hash== arg0.hash){return 0;}return (this.hash > arg0.hash) ? 1 : -1 ;
	}
}
