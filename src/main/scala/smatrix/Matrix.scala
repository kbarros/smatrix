package smatrix



// TODO: Specialize S#A for apply/update methods

abstract class Matrix[S <: Scalar : ScalarOps, +Repr[s <: Scalar] <: Matrix[s, Repr]]
    (val numRows: Int, val numCols: Int) extends MatrixDims { self: Repr[S] =>
  MatrixDims.checkDims(numRows, numCols)

  val scalar = implicitly[ScalarOps[S]]
  
  /** A description of the matrix type. */
  val description: String
  
  /** The matrix indices whose elements may be nonzero. For matrices with packed data buffers, the order
   * of these indices corresponds to the packing order.
   * TODO: Replace Iterable with Iterator
   */
  def definedIndices: Iterable[(Int, Int)]
  
  /** Gets the value of this matrix at given row and column indices.
   */
  def apply(i: Int, j: Int): S#A
  
  /** Sets the value of this matrix at given row and column indices.
   */
  def update(i: Int, j: Int, x: S#A)
    
  /** Assigns this matrix by applying a function to all matrix elements. For sparse matrices,
   * the transformation function must preserve zero. 
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

  /** Builds a new matrix by applying a function to all matrix elements. For sparse matrices,
   * the transformation must preserve zero.
   */
  // The parameter A2 is for type inference only
  def map[A2, S2 <: Scalar{type A=A2}, That[s <: Scalar] >: Repr[s] <: Matrix[s, That]]
      (f: S#A => S2#A)(implicit scalar2: ScalarOps[S2], mb: MatrixBuilder[S2, That]): That[S2] = {
    mb.map(this)(f)
  }

  /** This matrix scaled by a parameter.
   */
  def *[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.mul(_, x))
  }

  /** This matrix scaled by the inverse of a parameter.
   */
  def /[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.div(_, x))
  }

  /** The negation of this matrix.
   */
  def unary_-[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.neg(_))
  }

  /** The complex conjugate of this matrix.
   */
  def conj[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.conj(_))
  }

  /** The transpose of this matrix.
   */
  def tran[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    mb.transpose(this)
  }
  
  /** The hermitian conjugate (dagger) of this matrix.
   */
  def dag[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    tran(mb).transform(scalar.conj(_))
  }
  
  /** The dot product of this matrix with another, Tr(A B)
   */
  def dot[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2]]
      (that: M2[S])(implicit sb: ScalarBuilder[S], md: MatrixDotter[S, M1, M2]): S#A = {
    MatrixDims.checkDot(false, this, that)
    val accum = sb.build(1)
    md.dotTo(false, this, that, accum)
    val ret = scalar.read(accum, 0)
    accum.dispose()
    ret
  }
  
  /** The dot product of the conjugate of this matrix with another, Tr(A^dag B)
   */
  def dagDot[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2]]
      (that: M2[S])(implicit sb: ScalarBuilder[S], md: MatrixDotter[S, M1, M2]): S#A = {
    MatrixDims.checkDot(true, this, that)
    val accum = sb.build(1)
    md.dotTo(true, this, that, accum)
    val ret = scalar.read(accum, 0)
    accum.dispose()
    ret
  }

  /** Adds this matrix with another.
   */
  def +[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (that: M2[S])
        (implicit ma: MatrixPairAdder[S, M1, M2, M3],
                  mb: MatrixBuilder[S, M3]): M3[S] = {
    MatrixDims.checkAdd(this, that)
    val ret = mb.zeros(numRows, numCols)
    ma._1.addTo(scalar.one, this, ret)
    ma._2.addTo(scalar.one, that, ret)
    ret
  }

  /** Subtracts another matrix from this one.
   */
  def -[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (that: M2[S])
        (implicit ma: MatrixPairAdder[S, M1, M2, M3],
                  mb: MatrixBuilder[S, M3]): M3[S] = {
    MatrixDims.checkAdd(this, that)
    val ret = mb.zeros(numRows, numCols)
    ma._1.addTo(scalar.one, this, ret)
    ma._2.addTo(scalar.neg(scalar.one), that, ret)
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
  
  /**
   * Assigns `this := 0`.
   */
  def clear(): this.type = {
    transform(_ => scalar.zero)
  }
  
  /**
   * Assigns this matrix from another.
   */
  def :=[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2]]
        (that: M2[S])
        (implicit ma: MatrixAdder[S, M2, M1]): this.type = {
    MatrixDims.checkAssign(that, this)
    clear()
    ma.addTo(scalar.one, that, this)
    this
  }

  /** Assigns `this := alpha this`. 
   */
  def *=(x: S#A): this.type = {
    transform(scalar.mul(_, x))
  }

  /** Assigns `this := (1 / alpha) this`. 
   */
  def /=(x: S#A): this.type = {
    transform(scalar.div(_, x))
  }

  /** Assigns `this := this + that`.
   */
  def +=[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2]]
        (that: M2[S])
        (implicit ma: MatrixAdder[S, M2, M1]): this.type = {
    MatrixDims.checkAdd(this, that)
    ma.addTo(scalar.one, that, this)
    this
  }

  /** Assigns `this := this - that`.
   */
  def -=[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2]]
        (that: M2[S])
        (implicit ma: MatrixAdder[S, M2, M1]): this.type = {
    MatrixDims.checkAdd(this, that)
    ma.addTo(scalar.neg(scalar.one), that, this)
    this
  }

  /** Assigns `this := alpha that`.
   */
  def :=*[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2]]
        (alpha: S#A, that: M2[S])
        (implicit ma: MatrixAdder[S, M2, M1]): this.type = {
    MatrixDims.checkAssign(that, this)
    clear()
    ma.addTo(alpha, that, this)
    this
  }

  /** Assigns `this := this + alpha that`.
   */
  def +=*[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2]]
        (alpha: S#A, that: M2[S])
        (implicit ma: MatrixAdder[S, M2, M1]): this.type = {
    MatrixDims.checkAdd(that, this)
    ma.addTo(alpha, that, this)
    this
  }

  /** Assigns `this := m2 m3`.
   */
  def :=*[M1[s <: Scalar] >: Repr[s] <: Matrix[s, M1], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (m2: M2[S], m3: M3[S])
        (implicit mm: MatrixMultiplier[S, M2, M3, M1]): this.type = {
    gemm[M1, M2, M3](scalar.one, m2, m3, scalar.zero)
  }

  /** Assigns `this := this + m2 m3`.
   */
  def +=*[M1[s <: Scalar] >: Repr[s] <: Matrix[s, M1], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (m2: M2[S], m3: M3[S])
        (implicit mm: MatrixMultiplier[S, M2, M3, M1]): this.type = {
    gemm[M1, M2, M3](scalar.one, m2, m3, scalar.one)
  }

  /** Assigns `this := this - m2 m3`.
   */
  def -=*[M1[s <: Scalar] >: Repr[s] <: Matrix[s, M1], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (m2: M2[S], m3: M3[S])
        (implicit mm: MatrixMultiplier[S, M2, M3, M1]): this.type = {
    gemm[M1, M2, M3](scalar.neg(scalar.one), m2, m3, scalar.one)
  }

  /** Assigns `this := alpha m2 m3 + beta this`
   */
  def gemm[M1[s <: Scalar] >: Repr[s] <: Matrix[s, M1], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (alpha: S#A, m2: M2[S], m3: M3[S], beta: S#A)
        (implicit mm: MatrixMultiplier[S, M2, M3, M1]): this.type = {
    MatrixDims.checkMulTo(m2, m3, this)
    if (beta == scalar.zero)
      clear()
    else if (beta == scalar.one)
      ()
    else
      this *= beta
    mm.maddTo(alpha, m2, m3, this)
    this
  }

    
  /** The Frobenius norm squared of this matrix, Tr(A^dag A).
   */
  def norm2: S#A = {
    var ret = scalar.zero
    for ((i, j) <- definedIndices) {
      val x = this(i, j)
      ret = scalar.add(ret, scalar.mul(x, scalar.conj(x)))
    }
    ret
  }

  /** The sum of all matrix elements.
   */
  def sum: S#A = {
    var ret = scalar.zero
    for ((i, j) <- definedIndices) {
      ret = scalar.add(ret, this(i, j))
    }
    ret
  }
  
  /**
   * The matrix trace, \sum A_{ii}
   */
  def trace: S#A = {
    MatrixDims.checkTrace(this)
    var ret = scalar.zero
    for (i <- 0 until numRows) {
      ret = scalar.add(ret, this(i, i))
    }
    ret
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

abstract class MatrixAdder[S <: Scalar, M1[_ <: Scalar], M2[_ <: Scalar]] {
  /** Assigns `ret += alpha m`
   */
  def addTo(alpha: S#A, m: M1[S], ret: M2[S])
}

abstract class MatrixMultiplier[S <: Scalar, M1[_ <: Scalar], M2[_ <: Scalar], M3[s <: Scalar] <: Matrix[s, M3]] {
  /** Assigns `ret += alpha m1 m2`
   */
  def maddTo(alpha: S#A, m1: M1[S], m2: M2[S], ret: M3[S])
  
}

abstract class MatrixDotter[S <: Scalar, M1[_ <: Scalar], M2[_ <: Scalar]] {
  /** Assigns `ret = Tr(m1 m2) = \sum_{ij} m1_{ij} m2_{ji}`
   */
  def dotTo(dag1: Boolean, m1: M1[S], m2: M2[S], ret: RawData[S#Raw, S#Buf])
}

/** Used to resolve the return type for matrix addition and subtraction.  
 */
class MatrixPairAdder[S <: Scalar, M1[_ <: Scalar], M2[_ <: Scalar], M3[_ <: Scalar]]
  (implicit val _1: MatrixAdder[S, M1, M3], val _2: MatrixAdder[S, M2, M3])
