package smatrix


object Matrix {
  def addTo[S <: Scalar, M1[_] <: MatrixDims, M2[_] <: MatrixDims]
        (neg: Boolean, m: M1[S], ret: M2[S])
        (implicit ma: MatrixSingleAdder[S, M1, M2]) {
    MatrixDims.checkAddTo(m, ret)
    ma.addTo(neg, m, ret)
  }

  def maddTo[S <: Scalar, M1[_] <: MatrixDims, M2[_] <: MatrixDims, M3[_] <: MatrixDims]
        (m1: M1[S], m2: M2[S], ret: M3[S])
        (implicit mm: MatrixMultiplier[S, M1, M2, M3]) {
    MatrixDims.checkMulTo(m1, m2, ret)
    mm.maddTo(m1, m2, ret)
  }
}


// TODO: Specialize S#A for apply/update methods
// Eclipse crashes when "s" parameter is renamed, and file is saved 
trait Matrix[S <: Scalar, +Repr[s <: Scalar] <: Matrix[s, Repr]] extends MatrixDims { self: Repr[S] =>
  /** A description of the matrix type. */
  val description: String
  
  /** Scalar operations associated with the type of the matrix elements
   */
  val scalar: ScalarOps[S]
  
  /** Gets the value of this matrix at given row and column indices.
   */
  def apply(i: Int, j: Int): S#A
  
  /** Sets the value of this matrix at given row and column indices.
   */
  def update(i: Int, j: Int, x: S#A)
  
  /** Modifies this matrix by applying a function to each matrix element.
   */
  def transform(f: S#A => S#A): this.type
  
  /** Disposes the memory stored by this matrix, making this matrix invalid (optional).
   */
  def dispose() {}
  
  /** Creates a copy of this matrix.
   */
  def duplicate[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    mb.duplicate(this)
  }
  
  /** Converts this matrix to a dense matrix.
   */
  def toDense(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(numRows, numCols)
    for (j <- 0 until numCols; i <- 0 until numRows) { ret(i, j) = this(i, j) }
    ret
  }
  
  /** Builds a new matrix by applying a function to all elements of this matrix.
   */
  // The parameter A2 is for type inference only
  def map[A2, S2 <: Scalar{type A=A2}, That[s <: Scalar] >: Repr[s] <: Matrix[s, That]]
      (f: S#A => S2#A)(implicit scalar2: ScalarOps[S2], mb: MatrixBuilder[S2, That]): That[S2] = {
    mb.map(this)(f)
  }
  
  /** Returns this matrix scaled by a parameter.
   */
  def *[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.mul(_, x))
  }

  /** Returns this matrix scaled by the inverse of a parameter.
   */
  def /[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.div(_, x))
  }

  /** Returns the element-wise negation of this matrix.
   */
  def unary_-[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.neg(_))
  }

  /** Returns the element-wise conjugate of this matrix.
   */
  def conj[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.conj(_))
  }

  /** Returns the transpose of this matrix.
   */
  def tran[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    mb.transpose(this)
  }
  
  /** Returns the hermitian conjugate (dagger) of this matrix.
   */
  def dag[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    tran(mb).transform(scalar.conj(_))
  }
  
  /** Returns the Frobenius norm squared of this matrix (i.e., the square of the Euclidean norm of the vector of matrix elements).
   */
  def norm2: S#A = {
    var ret = scalar.zero
    for (i <- 0 until numRows; j <- 0 until numCols) {
      val x = this(i, j)
      ret = scalar.add(ret, scalar.mul(x, scalar.conj(x)))
    }
    ret
  }

  /** Sums all matrix elements.
   */
  def sum: S#A = {
    var ret = scalar.zero
    for (i <- 0 until numRows; j <- 0 until numCols) {
      ret = scalar.add(ret, this(i, j))
    }
    ret
  }
  
  /** Copies the matrix elements to an array in column major order (rows change most rapidly).
   */
  def toArray(implicit m: Manifest[S#A]): Array[S#A] = {
    val ret = new Array[S#A](numRows*numCols)
    for (j <- 0 until numCols;
         i <- 0 until numRows) {
      ret(i + j*numRows) = this(i, j)
    }
    ret
  }
  
  /** Takes the dot product of this row matrix and a column matrix.
   */
  def dot[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (that: M2[S])
        (implicit mm: MatrixMultiplier[S, M1, M2, M3],
                  mb: MatrixBuilder[S, M3]): S#A = {
    MatrixDims.checkDot(this, that)
    val m3 = mb.zeros(1, 1)
    mm.maddTo(this, that, m3)
    val ret = m3(0, 0)
    m3.dispose()
    ret
  }

  /** Multiplies this matrix with another.
   */
  def *[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (that: M2[S])
        (implicit mm: MatrixMultiplier[S, M1, M2, M3],
                  mb: MatrixBuilder[S, M3]): M3[S] = {
    MatrixDims.checkMul(this, that)
    val ret = mb.zeros(numRows, that.numCols)
    mm.maddTo(this, that, ret)
    ret
  }

  /** Adds this matrix with another.
   */
  def +[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (that: M2[S])
        (implicit ma: MatrixAdder[S, M1, M2, M3],
                  mb: MatrixBuilder[S, M3]): M3[S] = {
    MatrixDims.checkAdd(this, that)
    val ret = mb.zeros(numRows, numCols)
    ma.addTo(sub=false, this, that, ret)
    ret
  }

  /** Subtracts another matrix from this one.
   */
  def -[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (that: M2[S])
        (implicit ma: MatrixAdder[S, M1, M2, M3],
                  mb: MatrixBuilder[S, M3]): M3[S] = {
    MatrixDims.checkAdd(this, that)
    val ret = mb.zeros(numRows, numCols)
    ma.addTo(sub=true, this, that, ret)
    ret
  }
  
  /** Generates a string representation of this matrix.
   */
  override def toString = {
    val sb = new StringBuilder()
    val maxRows = 6
    val maxCols = 6
    val elemSpacing = 2
    
    var elemWidth = 8
    for (i <- 0 until math.min(numRows, maxRows)) {
      for (j <- 0 until math.min(numCols, maxCols)) {
        elemWidth = math.max(elemWidth, this(i, j).toString.size+elemSpacing)
      }
    }
    
    def writeStr(str: String) {
      val spaces = elemWidth-str.size
      val margin1 = spaces / 2
      val margin2 = spaces - margin1
      sb.append(Seq.fill(margin1)(' ').mkString)
      sb.append(str)
      sb.append(Seq.fill(margin2)(' ').mkString)
    }
    
    sb.append("[ %dx%d Matrix (%s) ]\n".format(numRows, numCols, description))
    for (i <- 0 until math.min(numRows, maxRows)) {
      for (j <- 0 until math.min(numCols, maxCols)) {
        writeStr(this(i, j).toString)
      }
      if (i == 0 && numCols > maxCols)
        sb.append(" ...")
      sb.append("\n")
    }
    if (numRows > maxRows) {
      writeStr(":")
    }
    sb.toString
  }
}


abstract class MatrixBuilder[S <: Scalar, Repr[_ <: Scalar]] {
  def zeros(numRows: Int, numCols: Int): Repr[S]
  def duplicate(m: Repr[S]): Repr[S]
  def transpose(m: Repr[S]): Repr[S]
  def map[S0 <: Scalar](m: Repr[S0])(f: S0#A => S#A): Repr[S]
}

abstract class MatrixSingleAdder[S <: Scalar, M1[_ <: Scalar], M2[_ <: Scalar]] {
  def addTo(neg: Boolean, m: M1[S], ret: M2[S])
}

class MatrixAdder[S <: Scalar, M1[_ <: Scalar], M2[_ <: Scalar], M3[_ <: Scalar]]
    (a13: MatrixSingleAdder[S, M1, M3], a23: MatrixSingleAdder[S, M2, M3]) {
  def addTo(sub: Boolean, m1: M1[S], m2: M2[S], ret: M3[S]) {
    a13.addTo(neg=false, m1, ret)
    a23.addTo(neg=sub, m2, ret)
  }
}

abstract class MatrixMultiplier[S <: Scalar, M1[_ <: Scalar], M2[_ <: Scalar], M3[_ <: Scalar]] {
  def maddTo(m1: M1[S], m2: M2[S], ret: M3[S])
}
