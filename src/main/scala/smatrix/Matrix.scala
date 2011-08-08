package smatrix


object MatrixDims {
  def checkDims(numRows: Int, numCols: Int) {
    require(numRows > 0 && numCols > 0,
        "Cannot build matrix with non-positive dimensions [%d, %d]".format(numRows, numCols))
  }
  
  def checkKey(m: MatrixDims, i: Int, j: Int) {
    require(0 <= i && i < m.numRows && 0 <= j && j < m.numCols,
        "Matrix indices out of bounds: [%d %d](%d %d)".format(m.numRows, m.numCols, i, j))
  }
  
  def checkAdd(m1: MatrixDims, m2: MatrixDims) {
    require(m1.numRows == m2.numRows && m1.numCols == m2.numCols,
        "Cannot add/subtract matrices of shape [%d, %d] +- [%d, %d]".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }
  def checkAddTo(m: MatrixDims, ret: MatrixDims) {
    require(m.numRows == ret.numRows && m.numCols == ret.numCols,
        "Cannot add/subtract matrices of shape +-[%d, %d] + [%d, %d] -> [%d, %d]".format(m.numRows, m.numCols, ret.numRows, ret.numCols, ret.numRows, ret.numCols))
  }

  def checkMul(m1: MatrixDims, m2: MatrixDims) {
    require(m1.numCols == m2.numRows,
        "Cannot multiply matrices of shape [%d, %d] * [%d, %d]".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }
  def checkMulTo(m1: MatrixDims, m2: MatrixDims, ret: MatrixDims) {
    checkMul(m1, m2)
    require(
        ret.numRows == m1.numRows &&
        m1.numCols == m2.numRows &&
        m2.numCols == ret.numCols, "Cannot multiply matrices of shape: [%d, %d] * [%d, %d] -> [%d, %d]".format(
            m1.numRows, m1.numCols, m2.numRows, m2.numCols, ret.numRows, ret.numCols))
  }
  
  def checkDot(m1: MatrixDims, m2: MatrixDims) {
    checkMul(m1, m2)
    require(m1.numRows == 1 && m2.numCols == 1,
        "Dot product expects row and column vectors, found [%d, %d] * [%d, %d]".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }
}


trait MatrixDims {
  def numRows: Int
  def numCols: Int
}


// TODO: Specialize S#A for apply/update methods

// Eclipse crashes when "s" parameter is renamed, and file is saved 

trait Matrix[S <: Scalar, +Repr[s <: Scalar] <: Matrix[s, Repr]] extends MatrixDims { self: Repr[S] =>
  val scalar: ScalarOps[S]
  
  def apply(i: Int, j: Int): S#A
  def update(i: Int, j: Int, x: S#A)
  def transform(f: S#A => S#A): this.type
  def dispose() {}
  
  def duplicate[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    mb.duplicate(this)
  }
  
  def toDense(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(numRows, numCols)
    for (j <- 0 until numCols; i <- 0 until numRows) { ret(i, j) = this(i, j) }
    ret
  }
  
  // The parameter A2 is for type inference only
  def map[A2, S2 <: Scalar{type A=A2}, That[s <: Scalar] >: Repr[s] <: Matrix[s, That]]
      (f: S#A => S2#A)(implicit scalar2: ScalarOps[S2], mb: MatrixBuilder[S2, That]): That[S2] = {
    mb.map(this)(f)
  }
  
  def *[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.mul(_, x))
  }

  def /[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.div(_, x))
  }

  def unary_-[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.neg(_))
  }

  def conj[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.conj(_))
  }

  def tran[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    mb.transpose(this)
  }
  
  def dag[That[s <: Scalar] >: Repr[s] <: Matrix[s, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    tran(mb).transform(scalar.conj(_))
  }

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

  def *[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (that: M2[S])
        (implicit mm: MatrixMultiplier[S, M1, M2, M3],
                  mb: MatrixBuilder[S, M3]): M3[S] = {
    MatrixDims.checkMul(this, that)
    val ret = mb.zeros(numRows, that.numCols)
    mm.maddTo(this, that, ret)
    ret
  }

  def +[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (that: M2[S])
        (implicit ma1: MatrixAdder[S, M1, M3],
                  ma2: MatrixAdder[S, M2, M3],
                  mb: MatrixBuilder[S, M3]): M3[S] = {
    MatrixDims.checkAdd(this, that)
    val ret = mb.zeros(numRows, numCols)
    ma1.addTo(neg=false, this, ret)
    ma2.addTo(neg=false, that, ret)
    ret
  }

  def -[M1[s <: Scalar] >: Repr[s], M2[s <: Scalar] <: Matrix[s, M2], M3[s <: Scalar] <: Matrix[s, M3]]
        (that: M2[S])
        (implicit ma1: MatrixAdder[S, M1, M3],
                  ma2: MatrixAdder[S, M2, M3],
                  mb: MatrixBuilder[S, M3]): M3[S] = {
    MatrixDims.checkAdd(this, that)
    val ret = mb.zeros(numRows, numCols)
    ma1.addTo(neg=false, this, ret)
    ma2.addTo(neg=true,  that, ret)
    ret
  }
  
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
    for (i <- 0 until math.min(numRows, maxRows)) {
      for (j <- 0 until math.min(numCols, maxCols)) {
        writeStr(this(i, j).toString)
      }
      if (i == 0 && numCols > maxCols)
        sb.append(" ... (%d Cols)".format(numCols))
      sb.append("\n")
    }
    if (numRows > maxRows) {
      writeStr(":")
      sb.append("\n")
      writeStr("(%d Rows)".format(numRows))
    }
    sb.toString
  }
}


trait MatrixBuilder[S <: Scalar, Repr[_ <: Scalar]] {
  def zeros(numRows: Int, numCols: Int): Repr[S]
  def duplicate(m: Repr[S]): Repr[S]
  def transpose(m: Repr[S]): Repr[S]
  def map[S0 <: Scalar](m: Repr[S0])(f: S0#A => S#A): Repr[S]
}

trait MatrixAdder[S <: Scalar, M1[_ <: Scalar], M2[_ <: Scalar]] {
  def addTo(neg: Boolean, m: M1[S], ret: M2[S])
}

trait MatrixMultiplier[S <: Scalar, M1[_ <: Scalar], M2[_ <: Scalar], M3[_ <: Scalar]] {
  def maddTo(m1: M1[S], m2: M2[S], ret: M3[S])
}
