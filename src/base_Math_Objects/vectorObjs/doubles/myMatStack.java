package base_Math_Objects.vectorObjs.doubles;
/**
 * Matrix stack structure, to mimic GL transformation stack
 * @author 7strb
 *
 */
public class myMatStack {
	public myMatrix[] s;
	public int top;
	 
	public myMatStack(int matStackMaxHeight){
		s = new myMatrix[matStackMaxHeight];
		for (int row = 0; row < matStackMaxHeight; ++row){s[row] = new myMatrix();}  	
		top = 0;        //point top of stack at index of base matrix
	}//stack constructor	
	
	public void initStackLocation(int idx){  this.initStackLocation(idx, true);}	
	public void initStackLocation(int idx, boolean toIdent){  s[idx].initMat(toIdent);}	
	/**
	 * add the current top of the matrix stack to the matrix stack in a higher position
	 */
	public void push(){ ++top; initStackLocation(top);  for (int row = 0; row < 4; ++row){ for (int col = 0; col < 4; ++col){  s[top].m[row][col] = s[top - 1].m[row][col]; }}}//push     	
	/**
	 * replace the current top of the matrix stack with a new matrix
	 * @param newTopMatrix
	 */
	public void replaceTop(myMatrix newTopMatrix){for (int row = 0; row < 4; ++row){ for (int col = 0; col < 4; ++col){ s[top].m[row][col] = newTopMatrix.m[row][col]; }}}//replaceTop	
	/**
	 * return the top of the matrix stack without popping
	 * @return
	 */
	public myMatrix peek(){ return s[top].clone();	}//peek	
	/**
	 * remove and return top matrix on stack
	 * @return
	 */
	public myMatrix pop(){
		myMatrix oldTop = new myMatrix();
		oldTop = s[top].clone();
		if (top > 0) {	top--;		initStackLocation(top+1,false);    }  //reinitialize stack			
		else {		System.out.println("stack pop error");}
		return oldTop;
	}	
	/**
	 * returns a string representation of this stack
	 */
	public String toString(){
		String result = "";
		for (int si = 0; si < top+1; si++){	result += "Stack[" + si + "] =\n[" + s[si].toString() + "]\n";}//for si
		return result;
	}//to String method
}//myMatStack
