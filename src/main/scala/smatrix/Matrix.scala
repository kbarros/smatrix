package smatrix



// TODO: Specialize S#A for apply/update methods
// Eclipse crashes when "s" parameter is renamed, and file is saved 
abstract class Matrix[S <: Scalar : ScalarOps, +Repr[s <: Scalar] <: Matrix[s, Repr]]
    (val numRows: Int, val numCols: Int) extends MatrixDims { self: Repr[S] =>
  MatrixDims.checkDims(numRows, numCols)

  val scalar = implicitly[ScalarOps[S]]
  
  /** A description of the matrix type. */
  val description: String
  
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
  
  /** The matrix indices whose elements may be nonzero
   */
  def definedIndices: Iterable[(Int, Int)] = for (i <- 0 until numRows; j <- 0 until numCols) yield (i, j)
  
  /** Creates a copy of this matrix.
   */
  def duplicate[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    mb.duplicate(this)
  }
  
  /** Converts this matrix to a dense matrix.
   */
  def toDense(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(numRows, numCols)
    for ((i, j) <- definedIndices) { ret(i, j) = this(i, j) }
    ret
  }

  /** Copies the matrix elements to an array in column major order (rows change most rapidly).
   */
  def toArray(implicit m: Manifest[S#A]): Array[S#A] = {
    val ret = Array.fill[S#A](numRows*numCols)(scalar.zero)
    for ((i, j) <- definedIndices) {
      ret(i + j*numRows) = this(i, j)
    }
    ret
  }

  /** Builds a new matrix by applying a function to all elements of this matrix.
   */
  // The parameter A2 is for type inference only
  def map[A2, S2 <: Scalar{type A=A2}, That[s <: Scalar] >: Repr[s] <: Matrix[s, That]]
      (f: S#A => S2#A)(implicit scalar2: ScalarOps[S2], mb: MatrixBuilder[S2, That]): That[S2] = {
    mb.map(this)(f)
  }
  
    
  /** Returns the Frobenius norm squared of this matrix (i.e., the square of the Euclidean norm of the vector of matrix elements).
   */
  def norm2: S#A = {
    var ret = scalar.zero
    for ((i, j) <- definedIndices) {
      val x = this(i, j)
      ret = scalar.add(ret, scalar.mul(x, scalar.conj(x)))
    }
    ret
  }

  /** Sums all matrix elements.
   */
  def sum: S#A = {
    var ret = scalar.zero
    for ((i, j) <- definedIndices) {
      ret = scalar.add(ret, this(i, j))
    }
    ret
  }
  
  /**
   * Returns the matrix trace, the sum of the diagonal elements
   */
  def trace: S#A = {
    MatrixDims.checkTrace(this)
    var ret = scalar.zero
    for (i <- 0 until numRows) {
      ret = scalar.add(ret, this(i, i))
    }
    ret
  }
  
  /**
   * Modifies this matrix by setting all elements to zero.
   */
  def clear(): this.type = {
    transform(_ => scalar.zero)
  }
  
  /** Returns this matrix scaled by a parameter.
   */
  def *[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.mul(_, x))
  }

  /** Modifies this matrix by scaling with a parameter. 
   */
  def *=(x: S#A): this.type = {
    transform(scalar.mul(_, x))
  }

  /** Returns this matrix scaled by the inverse of a parameter.
   */
  def /[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.div(_, x))
  }

  /** Modifies this matrix by scaling with a parameter inverse.
   */
  def /=(x: S#A): this.type = {
    transform(scalar.div(_, x))
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
  
  /** Takes the dot product of this row matrix and a column matrix.
   */
  def dot[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2]]
      (that: M2[S])
      (implicit sb: ScalarBuilder[S], md: MatrixDotter[S, M1, M2]): S#A = {
    MatrixDims.checkDot(this, that)
    val accum = sb.build(1)
    md.dotTo(this, that, accum)
    val ret = scalar.read(accum, 0)
    accum.dispose()
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
    mm.maddTo(scalar.one, this, that, ret)
    ret
  }

  /** Modifies this matrix to become the product of two parameter matrices.
   */
  def :=*[M1[s <: Scalar] >: Repr[s] <: Matrix[s, M1], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (m2: M2[S], m3: M3[S])
        (implicit mm: MatrixMultiplier[S, M2, M3, M1]): this.type = {
    MatrixDims.checkMulTo(m2, m3, this)
    this.clear()
    require((this ne m2) && (this ne m3), "Illegal aliasing of matrix product.")
    mm.maddTo(scalar.one, m2, m3, this)
    this
  }

  /** Accumulates the product of two parameter matrices into this one.
   */
  def +=*[M1[s <: Scalar] >: Repr[s] <: Matrix[s, M1], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (m2: M2[S], m3: M3[S])
        (implicit mm: MatrixMultiplier[S, M2, M3, M1]): this.type = {
    MatrixDims.checkMulTo(m2, m3, this)
    require((this ne m2) && (this ne m3), "Illegal aliasing of matrix product.")
    mm.maddTo(scalar.one, m2, m3, this)
    this
  }

  /** Modifies this matrix to become `alpha m2 m3 + beta this`
   */
  def gemm[M1[s <: Scalar] >: Repr[s] <: Matrix[s, M1], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (alpha: S#A, m2: M2[S], m3: M3[S], beta: S#A)
        (implicit mm: MatrixMultiplier[S, M2, M3, M1]): this.type = {
    MatrixDims.checkMulTo(m2, m3, this)
    require((this ne m2) && (this ne m3), "Illegal aliasing of matrix product.")
    mm.gemm(alpha, m2, m3, beta, this)
    this
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

  /** Accumulates another matrix into this one.
   */
  def +=[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2]]
        (that: M2[S])
        (implicit ma: MatrixSingleAdder[S, M2, M1]): this.type = {
    MatrixDims.checkAdd(this, that)
    ma.addTo(neg=false, that, this)
    this
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

  /** Accumulates the negative of another matrix into this one.
   */
  def -=[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2]]
        (that: M2[S])
        (implicit ma: MatrixSingleAdder[S, M2, M1]): this.type = {
    MatrixDims.checkAdd(this, that)
    ma.addTo(neg=true, that, this)
    this
  }

  /**
   * Modifies this matrix by copying from another one.
   */
  def :=[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2]]
        (that: M2[S])
        (implicit ma: MatrixSingleAdder[S, M2, M1]): this.type = {
    clear()
    this.+=[M1, M2](that)
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

abstract class MatrixMultiplier[S <: Scalar : ScalarOps, M1[_ <: Scalar], M2[_ <: Scalar], M3[s <: Scalar] <: Matrix[s, M3]] {
  val zero = implicitly[ScalarOps[S]].zero
  val one = implicitly[ScalarOps[S]].one
  
  /** ret := alpha m1 m2 + beta ret
   */
  def gemm(alpha: S#A, m1: M1[S], m2: M2[S], beta: S#A, ret: M3[S]) {
    if (beta == zero)
      ret.clear()
    else if (beta == one)
      ()
    else {
      ret *= beta
    }
    maddTo(alpha, m1, m2, ret)
  }
  
  /** ret += alpha m1 m2
   */
  def maddTo(alpha: S#A, m1: M1[S], m2: M2[S], ret: M3[S])
  
}

abstract class MatrixDotter[S <: Scalar, M1[_ <: Scalar], M2[_ <: Scalar]] {
  def dotTo(m1: M1[S], m2: M2[S], ret: RawData[S#Raw, S#Buf])
}
