package base_Math_Objects.matrixObjs.doubles;

import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.doubles.myVector;

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
	public void push(){ 
		++top; 
		initStackLocation(top);  
		for (int row = 0; row < s[top].m.length; ++row){ 
			for (int col = 0; col < s[top].m[row].length; ++col){  
				s[top].m[row][col] = s[top - 1].m[row][col]; 
			}
		}
	}//push     	
	/**
	 * replace the current top of the matrix stack with a new matrix
	 * @param newTopMatrix
	 */
	public void replaceTop(myMatrix newTopMatrix){
		for (int row = 0; row < s[top].m.length; ++row){ 
			for (int col = 0; col < s[top].m[row].length; ++col){  
				s[top].m[row][col] = newTopMatrix.m[row][col]; 
			}
		}
	}//replaceTop	
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
	 * Build a translation transformation matrix and multiply it against the top matrix in the stack
	 * @param tx
	 * @param ty
	 * @param tz
	 */
	public void translate(double tx, double ty, double tz) { 
		//build and push onto stack the translation matrix
		myMatrix TransMat = new myMatrix();
		//set the 4th column vals to be the translation coordinates
		TransMat.setValByIdx(0,3,tx);
		TransMat.setValByIdx(1,3,ty);
		TransMat.setValByIdx(2,3,tz);
		updateCTM(TransMat);
	}//translate method
	
	/**
	 * Build a scale transformation matrix and multiply it against the top matrix in the stack
	 * @param sx
	 * @param sy
	 * @param sz
	 */
	public void scale(double sx, double sy, double sz) {
		//build and push onto stack the scale matrix
		myMatrix ScaleMat = new myMatrix();
		//set the diagonal vals to be the scale coordinates
		ScaleMat.setValByIdx(0,0,sx);
		ScaleMat.setValByIdx(1,1,sy);
		ScaleMat.setValByIdx(2,2,sz);
		updateCTM(ScaleMat);
	}//scale method

	/**
	*  Builds a rotation matrix to be in "angle" degrees CCW around the axis given by ax,ay,az
	*  and multiples this matrix against the CTM
	*/
	public void rotate(double angle, double ax, double ay, double az) { 
		// build and add to top of stack the rotation matrix
		double angleRad = angle * MyMathUtils.DEG_TO_RAD;
		myMatrix RotMat = new myMatrix();
		myMatrix RotMatrix1 = new myMatrix();      //translates given axis to x axis
		myMatrix RotMatrix2 = new myMatrix();      //rotation around x axis by given angle
		myMatrix RotMatrix1Trans = new myMatrix();
	  
		myVector axisVect, axisVectNorm, bVect, bVectNorm, cVect, cVectNorm, normVect;
		//first build rotation matrix to rotate ax,ay,az to lie in line with x axis		
		axisVect = new myVector(ax,ay,az);
		axisVectNorm = axisVect._normalized();
	  
		if (ax == 0) { 	normVect = new myVector(1,0,0);} 
		else {			normVect = new myVector(0,1,0);}
		bVect = axisVectNorm._cross(normVect);
		bVectNorm = bVect._normalized();
	  
		cVect = axisVectNorm._cross(bVectNorm);
		cVectNorm = cVect._normalized();
	  
		RotMatrix1.setValByRow(0,axisVectNorm);
		RotMatrix1.setValByRow(1,bVectNorm);
		RotMatrix1.setValByRow(2,cVectNorm);
		
		RotMatrix1Trans = RotMatrix1.transpose();
		//second build rotation matrix to rotate around x axis by angle
		//need to set 1,1 ; 1,2 ; 2,1 ; and 2,2 to cos thet, neg sine thet, sine thet, cos thet, respectively
		double cosVal = (Math.cos(angleRad)), sinVal =(Math.sin(angleRad)); 
		RotMatrix2.setValByIdx(1,1,cosVal);
		RotMatrix2.setValByIdx(1,2,-sinVal);
		RotMatrix2.setValByIdx(2,1,sinVal);
		RotMatrix2.setValByIdx(2,2,cosVal);
		//lastly, calculate full rotation matrix

		myMatrix tmp = RotMatrix2.multMat(RotMatrix1);
		RotMat = RotMatrix1Trans.multMat(tmp);
		updateCTM(RotMat);
	}//rotate method
	
	
	/**
	 * Update the current top matrix
	 * @param _mat
	 */
	private void updateCTM(myMatrix _mat){		
		myMatrix CTM = peek();
		replaceTop(CTM.multMat(_mat));
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
