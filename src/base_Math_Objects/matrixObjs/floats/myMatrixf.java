package base_Math_Objects.matrixObjs.floats;

import base_Math_Objects.vectorObjs.floats.myVectorf;

public class myMatrixf {
	public float[][] m;
	
	public myMatrixf(){  m = new float[4][4]; initMat();}
	public myMatrixf(float[][] _m) {
		m = new float[_m.length][_m[0].length];
		for(int row=0;row<m.length;++row) {
			for(int col=0;col<m[row].length;++col) {m[row][col]=_m[row][col];}	
		}
	}
	public void initMat(){  this.initMat(true);}	
	/**
	 * initialize this matrix to be identity matrix or all zeros
	 * @param toIdentity whether should be identity or not
	 */
	public void initMat(boolean toIdentity){
		for (int row = 0; row < m.length; ++row){for (int col = 0; col < m[row].length; ++col){m[row][col] = 0;}}
		if (toIdentity) {for(int diag = 0;diag<m.length;++diag) {m[diag][diag] = 1.0f;}}
	} 
	
	/**
	 * Returns identity matrix
	 * @return
	 */
	public static final myMatrixf Identity() {return new myMatrixf();}
	
	/**
	 * multiplies this matrix by b, in order: [this] x [b] returns result in result
	 * @param b
	 * @return
	 */
	public myMatrixf multMat(myMatrixf b){
		float resultVal = 0;
		myMatrixf result = new myMatrixf();
		for (int row = 0; row < this.m.length; ++row){
			for (int col = 0; col < this.m[row].length; ++col){
				for (int k = 0; k < b.m.length; k++){
					resultVal += this.m[row][k] * b.getValByIdx(k,col);
				}
				result.setValByIdx(row,col,resultVal); 
				resultVal = 0;
			}
		}
		return result;  
	}//mult method
		
	/**
	 * multiplies this matrix by vertex b, in order: [this] x [b]
	 * returns result vertex in result
	 * @param b
	 * @return
	 */
	public float[] multVert(float[] b){
		float resultVal;
		float[] result = new float[]{0,0,0,0};
		for (int row = 0; row < this.m.length; ++row){
			resultVal = 0;
			for (int col = 0; col < this.m[row].length; ++col){resultVal += this.m[row][col] * b[col];}	
			result[row] = resultVal;
		}//for row
		return result;  
	}//mult method
	
	
	/**
	 * Multiply ever element of this matrix by passed scale amount
	 * @param scl
	 */
	public void scaleMatrix(float scl) {
		for (int row = 0; row < this.m.length; ++row){
			for (int col = 0; col < this.m[row].length; ++col){ this.m[row][col] *= scl;}
		}
	}//scaleMatrix
	

	/**
	 * Finds the inverse of this matrix, if it is invertible. Otherwise, returns Identity matrix
	 * @return
	 */
	public myMatrixf inverse(){
		//first find determinant.  If is 0 then no inverse
		float det = determinant();
		if(det == 0) {
			System.out.println("Uninvertible matrix -> determinant == 0. Aborting");
			return new myMatrixf();
		}
		//Start with adjoint
		myMatrixf result = adjoint();
		//now multiply by 1/det
		result.scaleMatrix(1.0f/det);
		return result; 
	}//method invert
	
	/**
	 * Return the adjoint of this matrix - the transpose of the cofactor matrix
	 * @return
	 */
	public myMatrixf adjoint(){
		//transpose of cofactor matrix
		float[][] res = myMatrixf.AdjointMatrix(this.m);
		myMatrixf newMat = new myMatrixf(res);
		return newMat;
	}
	
	/**
	 * Return the Ajoint of the passed matrix of floats - the transpose of the cofactor matrix
	 * @param M
	 * @return
	 */
	public synchronized static float[][] AdjointMatrix(float[][] M){
		float[][] coFactor = myMatrixf.CoFactorMatrix(M);
		float[][] res = myMatrixf.TransposeMatrix(coFactor);
		return res;	
	}
	
	/**
	 * returns the transpose of this matrix - also inverse if rotation matrix
	 * @return
	 */
	public myMatrixf transpose(){ 
		float[][] res = myMatrixf.TransposeMatrix(this.m);
		myMatrixf newMat = new myMatrixf(res);
		return newMat; 
	}//transpose method
	
	/**
	 * returns the transpose of the passed float matrix
	 * @return
	 */
	public synchronized static float[][] TransposeMatrix(float[][] M){
		float[][] res = new float[M.length][M[0].length];
		for(int row=0;row<M.length;++row) {
			for(int col=0;col<M[0].length;++col) {
				res[col][row] = M[row][col];
			}
		}		
		return res;
	}//TransposeMatrix
	
	/**
	 * Return the cofactor matrix of this matrix
	 * @return
	 */
	public myMatrixf coFactorMatrix() {
		float[][] res = myMatrixf.CoFactorMatrix(this.m);
		myMatrixf newMat = new myMatrixf(res);
		return newMat;
	}
	
	/**
	 * Returns the cofactor matrix of the passed float matrix
	 * @param M
	 * @return
	 */
	public synchronized static float[][] CoFactorMatrix(float[][] M){
		float[][] res = new float[M.length][M[0].length];
		float[][] scratch = new float[M.length-1][M[0].length-1];
		float s = 1.0f, colS;
		for(int col=0;col<res[0].length;++col) {
			colS = s;
			for(int row=0;row<res.length;++row) {
		         int sRow = 0;
		         for (int mRow=0;mRow<res.length;++mRow) {
		            if (mRow == row) {continue;}
		            int sCol = 0;
		            for (int mCol=0;mCol<res[row].length;++mCol) {
		               if (mCol == col){continue;}
		               scratch[sRow][sCol] = M[mRow][mCol];
		               ++sCol;
		            }
		            ++sRow;
		         }		         
		         //adding to prevent negative 0
		         res[row][col] += s * myMatrixf.Determinant(scratch);				
		         s *= -1.0f;
			}		
			//flip the starting s for next column
			s = colS * -1.0f;
		}			
		return res;
	}//CoFactorMatrix	
	
	/**
	 * calculates the determinant of this 4x4 matrix
	 * @return
	 */
	public float determinant() {
		return myMatrixf.Determinant(this.m);
	}

	/**
	 * calculates the determinant of a Matrix of floats
	 * @param M n x n matrix - don't over do it
	 * @param coFactors n x n matrix (OUT) resultant coFactors
	 * @return determinant
	 */
	public synchronized static float Determinant(float[][] M){ 
		float sum=0, s; 
		if (M.length == 2) {//2 x 2 mat
			return M[0][0] * M[1][1] - M[1][0] * M[0][1];
		}
		for(int col=0;col < M.length;++col){ 	
			//construct minor w/respect to elements of top row of M
			float[][] minor= new float[M.length-1][M.length-1];
			for(int mCol=0;mCol<M.length;++mCol){
				if(mCol==col) {continue;}
				int mColIdx = (mCol<col)? mCol : mCol-1;
				for(int mRow=1;mRow<M.length;++mRow){
					minor[mRow-1][mColIdx] = M[mRow][mCol];
				}
			}	
			s = (col%2==0) ? 1.0f : -1.0f;
			sum += s * M[0][col] * Determinant(minor); 										
		}
		return(sum); //returns determinant value. once stack is finished, returns final determinant.
	}//Determinant
	
	
	/**
	 * Return this matrix as an array in row or column major order based on passed boolean
	 * @return
	 */
	public float[] getMatAsArray(boolean isRowMajor) {
		int numElemsPerRow = this.m.length;
		int numElemsPerCol = this.m[0].length;
		float[] res = new float[numElemsPerRow*numElemsPerCol];
		
		if (isRowMajor) {
			//contiguous across each row
			for(int row=0; row<numElemsPerRow;++row) {
				int rowIncr = row * numElemsPerCol;
				for (int col=0;col< numElemsPerCol;++col) {
					res[rowIncr + col] = this.m[row][col];	
				}
			}			
		} else {
			//colMajor : contiguous across each column
			for (int col=0;col< numElemsPerCol;++col) {
				int colIncr = col * numElemsPerRow;
				for(int row=0; row<numElemsPerRow;++row) {
					res[row + colIncr] = this.m[row][col];	
				}
			}			
		}		
		return res;
	}//getMatAsArray
		
	public float getValByIdx(int row, int col){   return m[row][col]; }  
	public void setValByIdx(int row, int col, float val){	   m[row][col] = val; }
	
		//writes over the first 3 cols of row given by row var with vector vals  
	public void setValByRow(int row, myVectorf vect){ this.m[row][0] = vect.x; this.m[row][1] = vect.y;      this.m[row][2] = vect.z; }//setValByRow
	 
	//makes a deep copy of this matrix, which it returns
	public myMatrixf clone(){
		myMatrixf newMat = new myMatrixf();
		for (int row = 0; row < this.m.length; ++row){ for (int col = 0; col < this.m[row].length; ++col){newMat.m[row][col] = this.m[row][col];}}
		return newMat;
	}
	
	public String toString(){
	    String result = "", tmp2str = "",tmpString;
	    for (int row = 0; row < this.m.length; ++row){
	    	result += "[";
	    	for (int col = 0; col < this.m[row].length; ++col){   tmp2str = "" + m[row][col]; if (col != this.m[row].length-1) {tmp2str += ", ";} result += tmp2str;}
	    	tmpString = "]";  if (row != this.m.length-1) { tmpString += "\n"; }
	    	result += tmpString;
    	}//for row	   
	    return result;
	}//toString method
  
}//class matrixf