package base_Math_Objects.vectorObjs.doubles;

public class myMatrix {	  
	public double[][] m;
	
	public myMatrix(){  m = new double[4][4]; initMat();}
	public myMatrix(double[][] _m) {
		m = new double[4][4];
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
		if (toIdentity) {for(int diag = 0;diag<m.length;++diag) {m[diag][diag] = 1.0;}}
	} 
	
	/**
	 * multiplies this matrix by b, in order: [this] x [b] returns result in result
	 * @param b
	 * @return
	 */
	public myMatrix multMat(myMatrix b){
		double resultVal = 0;
		myMatrix result = new myMatrix();
		for (int row = 0; row < this.m.length; ++row){
			for (int col = 0; col < this.m[row].length; ++col){
				for (int k = 0; k < this.m[row].length; k++){
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
	public double[] multVert(double[] b){
		double resultVal;
		double[] result = new double[]{0,0,0,0};
		for (int row = 0; row < 4; ++row){
			resultVal = 0;
			for (int col = 0; col < 4; ++col){resultVal += this.m[row][col] * b[col];}	
			result[row] = resultVal;
		}//for row
		return result;  
	}//mult method

	
	/**
	 * Multiply ever element of this matrix by passed scale amount
	 * @param scl
	 */
	public void scaleMatrix(double scl) {
		for (int row = 0; row < this.m.length; ++row){
			for (int col = 0; col < this.m[row].length; ++col){ this.m[row][col] *= scl;}
		}
	}//scaleMatrix
	
	/**
	 * Finds the inverse of this matrix, if it is invertible. Otherwise, returns Identity matrix
	 * @return
	 */
	public myMatrix inverse(){
		//first find determinant.  If is 0 then no inverse
		double det = determinant();
		if(det == 0) {
			System.out.println("Uninvertible matrix -> determinant == 0. Aborting");
			return new myMatrix();
		}
		//Start with adjoint
		myMatrix result = adjoint();
		//now multiply by 1/det
		result.scaleMatrix(1.0/det);
		return result;    
	}//method invert
	
	/**
	 * Finds the inverse of this matrix, if it is invertible. Otherwise, returns Identity matrix
	 * @return
	 */
	public myMatrix inverseOld(){
		myMatrix result = this.InvertMe();
		return result;    
	}//method invert
	/**
	 * Find the inverse of the passed matrix, by finding the determinant, checking value to see if invertible, 
	 * and then finding the Adjoint, and scaling it by 1/det. If not invertible, returns M sized Identity
	 * @param M
	 * @return
	 */
	public synchronized static double[][] InvertMatrix(double[][] M){		
		double det = myMatrix.Determinant(M);
		if(det == 0) {
			System.out.println("Uninvertible matrix -> determinant == 0. Aborting");
			double[][] res = new double[M.length][M[0].length];
			for(int i=0;i<M.length;++i) {	res[i][i] = 1.0;}
			return res;
		}
		double[][] res = myMatrix.AdjointMatrix(M);
		//scale by determinant
		for(int row=0;row<res.length;++row) {
			for(int col=0;col<res[row].length;++col) {res[row][col] /= det;}
		}		
		return res;
	}//InvertMatrix
	
	/**
	 * Return the adjoint of this matrix - the transpose of the cofactor matrix
	 * @return
	 */
	public myMatrix adjoint(){
		//transpose of cofactor matrix
		double[][] res = myMatrix.AdjointMatrix(this.m);
		myMatrix newMat = new myMatrix(res);
		return newMat;
	}
	
	/**
	 * Return the Ajoint of the passed matrix of doubles - the transpose of the cofactor matrix
	 * @param M
	 * @return
	 */
	public synchronized static double[][] AdjointMatrix(double[][] M){
		double[][] coFactor = myMatrix.CoFactorMatrix(M);
		double[][] res = myMatrix.TransposeMatrix(coFactor);
		return res;	
	}
	
	
	/**
	 * returns the transpose of this matrix - also inverse if rotation matrix
	 * @return
	 */
	public myMatrix transpose(){ 
		double[][] res = myMatrix.TransposeMatrix(this.m);
		myMatrix newMat = new myMatrix(res);
		return newMat; 
	}//transpose method
	
	/**
	 * returns the transpose of the passed double matrix
	 * @return
	 */
	public synchronized static double[][] TransposeMatrix(double[][] M){
		double[][] res = new double[M.length][M[0].length];
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
	public myMatrix coFactorMatrix() {
		double[][] res = myMatrix.CoFactorMatrix(this.m);
		myMatrix newMat = new myMatrix(res);
		return newMat;
	}
	
	/**
	 * Returns the cofactor matrix of the passed float matrix - find determinant of every minor
	 * @param M
	 * @return
	 */
	public synchronized static double[][] CoFactorMatrix(double[][] M){
		double[][] res = new double[M.length][M[0].length];
		double[][] scratch = new double[M.length-1][M[0].length-1];
		double s = 1.0, colS;		
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
		         res[row][col] += s * myMatrix.Determinant(scratch);				
		         s *= -1.0;
			}
			//flip the starting s for next column
			s = colS * -1.0;
		}		
		return res;
	}//CoFactorMatrix
	
	
	/**
	 * calculates the determinant of this 4x4 matrix
	 * @return
	 */
	public double determinant() {
		return myMatrix.Determinant(this.m);
	}
	
	/**
	 * calculates the determinant of a Matrix of doubles
	 * @param M n x n matrix - don't over do it
	 * @return
	 */

	public synchronized static double Determinant(double[][] M){ 
		double sum=0, s;
		if (M.length ==2) {//2 x 2 mat
			return M[0][0] * M[1][1] - M[1][0] * M[0][1];
		}
		for(int col=0;col < M.length;++col){ 	
			//construct minor w/respect to elements of top row of M
			double[][] minor= new double[M.length-1][M.length-1];
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
	public double[] getMatAsArray(boolean isRowMajor) {
		int numElemsPerRow = this.m.length;
		int numElemsPerCol = this.m[0].length;
		double[] res = new double[numElemsPerRow*numElemsPerCol];
		
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
	
	
	//------------- inversion code
	/**
	 * Invert this matrix - replaced by above determinant + cofactor
	 * @return
	 */
	private myMatrix InvertMe(){
		double[] tmp = new double[12];   // temp array for pairs
		double[] src = new double[16];   // array of transpose source matrix
		double[] dst = new double[16];   //destination matrix, in array form
		double det = 0;       // determinant 
		myMatrix dstMat = new myMatrix();
		//convert this matrix to array form and set up source vector
		for(int row = 0; row < this.m.length; ++row){ 
			for(int col = 0; col < this.m[row].length; ++col){ 
				src[(4*col) + row] = this.m[row][col];	
			}
		}
		
		// calculate pairs for first 8 elements (cofactors)
		tmp[0] = src[10] * src[15];
		tmp[1] = src[11] * src[14];
		tmp[2] = src[9] * src[15];
		tmp[3] = src[11] * src[13];
		tmp[4] = src[9] * src[14];
		tmp[5] = src[10] * src[13];
		tmp[6] = src[8] * src[15];
		tmp[7] = src[11] * src[12];
		tmp[8] = src[8] * src[14];
		tmp[9] = src[10] * src[12];
		tmp[10] = src[8] * src[13];
		tmp[11] = src[9] * src[12];
		
		// calculate first 8 elements (cofactors) 	  
		dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
		dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
		dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
		dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
		dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
		dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
		dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
		dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
		dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
		dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
		dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
		dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
		dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
		dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
		dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
		dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
		
		// calculate pairs for second 8 elements (cofactors) 
		
		tmp[0] = src[2]*src[7];
		tmp[1] = src[3]*src[6];
		tmp[2] = src[1]*src[7];
		tmp[3] = src[3]*src[5];
		tmp[4] = src[1]*src[6];
		tmp[5] = src[2]*src[5];
		tmp[6] = src[0]*src[7];
		tmp[7] = src[3]*src[4];
		tmp[8] = src[0]*src[6];
		tmp[9] = src[2]*src[4];
		tmp[10] = src[0]*src[5];
		tmp[11] = src[1]*src[4];
		
		// calculate second 8 elements (cofactors)
		
		dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
		dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
		dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
		dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
		dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
		dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
		dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
		dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
		dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
		dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
		dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
		dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
		dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
		dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
		dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
		dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
		
		// calculate determinant 		
		det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];
		if(Math.abs(det) > .0000001){		
			for (int j = 0; j < dst.length; ++j){ dst[j] /= det; }	    
			 //convert dst array to matrix
			for(int row = 0; row < dstMat.m.length; ++row){ 
				for(int col = 0; col < dstMat.m[row].length; ++col){  
					dstMat.m[row][col] = dst[(dstMat.m[row].length*row) + col];
				}
			}
		} else {
			System.out.println("uninvertible matrix -> det == 0");
		}
		return dstMat;
	}//invertme code	
	//end------------inversion code
	
	public double getValByIdx(int row, int col){   return m[row][col]; }  
	public void setValByIdx(int row, int col, double val){	   m[row][col] = val; }
	
	/**
	 * writes over the first 3 cols of row given by row var with vector vals  
	 * @param row
	 * @param vect
	 */
	public void setValByRow(int row, myVector vect){ this.m[row][0] = vect.x; this.m[row][1] = vect.y;      this.m[row][2] = vect.z; }//setValByRow
	 
	/**
	 * makes a deep copy of this matrix, which it returns
	 */
	public myMatrix clone(){
		myMatrix newMat = new myMatrix();
		for (int row = 0; row < this.m.length; ++row){ for (int col = 0; col < this.m[row].length; ++col){newMat.m[row][col] = this.m[row][col];}}
		return newMat;
	}
	
	public String toString(){
	    String result = "", tmp2str = "",tmpString;
	    for (int row = 0; row < 4; ++row){
	    	result += "[";
	    	for (int col = 0; col < this.m[row].length; ++col){   tmp2str = "" + m[row][col]; if (col != this.m[row].length-1) {tmp2str += ", ";} result += tmp2str;}
	    	tmpString = "]";  if (row != this.m.length-1) { tmpString += "\n"; }
	    	result += tmpString;
    	}//for row	   
	    return result;
	}//toString method
  
}//class matrix