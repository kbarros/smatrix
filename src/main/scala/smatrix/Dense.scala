package smatrix


trait DenseSlice
object :: extends DenseSlice


object Dense extends DenseBuilders with DenseAdders with DenseMultipliers with DenseLapackImplicits


class Dense[S <: Scalar : ScalarOps : ScalarBuilder : Netlib](numRows: Int, numCols: Int)
    extends Matrix[S, Dense](numRows, numCols) {
  override val description = "Dense"
  val netlib: Netlib[S] = implicitly[Netlib[S]]
  val data: RawData[S#Raw, S#Buf] = implicitly[ScalarBuilder[S]].build(numRows*numCols)
  
  def index(i: Int, j: Int) = {
    MatrixDims.checkKey(this, i, j)
    i + j*numRows // fortran column major convention
  }
  
  override def apply(i: Int, j: Int): S#A = scalar.read(data, index(i, j))
  def apply(i: Int): S#A = {
    if (numRows == 1)
      this(0, i)
    else if (numCols == 1)
      this(i, 0)
    else {
      require(numRows == 1 || numCols == 1, "Cannot apply a single index to non-vector matrix of shape [%d, %d].".format(numRows, numCols))
      sys.error("")
    }
  }
  def apply(i: Int, _slice: DenseSlice)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(1, numCols)
    for (j <- 0 until numCols) ret(0, j) = this(i, j)
    ret
  }
  def apply(_slice: DenseSlice, j: Int)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(numRows, 1)
    for (i <- 0 until numRows) ret(i, 0) = this(i, j)
    ret
  }
  
  override def update(i: Int, j: Int, x: S#A): Unit = scalar.write(data, index(i, j), x)
  def update(i: Int, x: S#A): Unit = {
    if (numRows == 1)
      this(0, i) = x
    else if (numCols == 1)
      this(i, 0) = x
    else {
      require(numRows == 1 || numCols == 1, "Cannot apply a single index to non-vector matrix of shape [%d, %d].".format(numRows, numCols))
      sys.error("")
    }
  }
  def update[That[s <: Scalar] <: Matrix[s, That]](i: Int, _slice: DenseSlice, that: That[S]) {
    require(that.numRows == 1 && numCols == that.numCols, "Cannot perform matrix assignment of shape: [%d, %d] -> [%d, %d](%d, ::).".format(
      that.numRows, that.numCols, numRows, numCols, i))
    for (j <- 0 until numCols) this(i, j) = that(0, j)
  }
  def update[That[s <: Scalar] <: Matrix[s, That]](_slice: DenseSlice, j: Int, that: That[S]) {
    require(that.numCols == 1 && numRows == that.numRows, "Cannot perform matrix assignment of shape: [%d, %d] -> [%d, %d](::, %d).".format(
      that.numRows, that.numCols, numRows, numCols, j))
    for (i <- 0 until numRows) this(i, j) = that(i, 0)
  }
   
  override def transform(f: S#A => S#A): this.type = {
    for (i <- 0 until numRows*numCols) { scalar.write(data, i, f(scalar.read(data, i))) }
    this
  }
}


trait DenseBuilders {
  implicit def denseBuilder[S <: Scalar : ScalarOps : ScalarBuilder : Netlib] = new MatrixBuilder[S, Dense] {
    def zeros(numRows: Int, numCols: Int) = {
      new Dense[S](numRows, numCols)
    }
        
    def duplicate(m: Dense[S]): Dense[S] = {
      val ret = zeros(m.numRows, m.numCols)
      m.data.copyTo(ret.data)
      ret
    }

    def transpose(m: Dense[S]): Dense[S] = {
      val ret = zeros(m.numCols, m.numRows)
      for (i <- 0 until ret.numRows; j <- 0 until ret.numCols) { ret(i, j) = m(j, i) } 
      ret
    }
    
    def map[S0 <: Scalar](m: Dense[S0])(f: S0#A => S#A): Dense[S] = {
      val ret = zeros(m.numRows, m.numCols)
      for (i <- 0 until ret.numRows; j <- 0 until ret.numCols) { ret(i, j) = f(m(i, j)) } 
      ret
    }
  }
}


trait DenseAdders {
  implicit def denseDenseAdder[S <: Scalar] = new MatrixAdder[S, Dense, Dense] {
    def addTo(alpha: S#A, m: Dense[S], ret: Dense[S]) = {
      MatrixDims.checkAdd(m, ret)
      for (i <- 0 until ret.numRows;
           j <- 0 until ret.numCols) {
        ret.scalar.maddTo(m.data, m.index(i, j), alpha, ret.data, ret.index(i, j))
      }
    }
  }
  
  implicit def denseDensePairAdder[S <: Scalar] = new MatrixPairAdder[S, Dense, Dense, Dense]
}


trait DenseMultipliers {
  implicit def denseDenseMultiplier[S <: Scalar] = new MatrixMultiplier[S, Dense, Dense, Dense] {
    def maddTo(alpha: S#A, m1: Dense[S], m2: Dense[S], ret: Dense[S]) {
      MatrixDims.checkMulTo(m1, m2, ret)
      val s = ret.scalar
      if (alpha == s.zero) {
        // no-op
      }
      else if (ret.netlib == null || (ret.numCols min ret.numRows) <= 4) {
        if (alpha == s.one) {
          for (j <- 0 until ret.numCols;
               k <- 0 until m1.numCols;
               i <- 0 until ret.numRows) {
            s.maddTo(m1.data, i+k*m1.numRows, m2.data, k+j*m2.numRows, ret.data, i+j*ret.numRows)
          }
        }
        else if (ret.numRows > ret.numCols) {
          for (j <- 0 until ret.numCols;
               k <- 0 until m1.numCols;
               val alpha_m2_kj = s.mul(alpha, s.read(m2.data, k+j*m2.numRows)); 
               i <- 0 until ret.numRows) {
            s.maddTo(m1.data, i+k*m1.numRows, alpha_m2_kj, ret.data, i+j*ret.numRows)
          }
        }
        else {
          for (i <- 0 until ret.numRows;
               k <- 0 until m1.numCols;
               val alpha_m1_ik = s.mul(alpha, s.read(m1.data, i+k*m1.numRows)); 
               j <- 0 until ret.numCols) {
            s.maddTo(m2.data, k+j*m2.numRows, alpha_m1_ik, ret.data, i+j*ret.numRows)
          }
        }
      }
      else {
        ret.netlib.gemm(Netlib.CblasColMajor, Netlib.CblasNoTrans, Netlib.CblasNoTrans,
            ret.numRows, ret.numCols, // dimension of return matrix
            m1.numCols, // dimension of summation index
            alpha, // alpha
            m1.data.buffer, m1.numRows, // A matrix
            m2.data.buffer, m2.numRows, // B matrix
            s.zero, // beta
            ret.data.buffer, ret.numRows // C matrix
        )
      }
    }
  }
  
  implicit def denseDenseDotter[S <: Scalar] = new MatrixDotter[S, Dense, Dense] {
    def dotTo(m1: Dense[S], m2: Dense[S], ret: RawData[S#Raw, S#Buf]) {
      MatrixDims.checkDot(m1, m2)
      val s = m1.scalar
      s.write(ret, 0, s.zero)
      var i = 0
      while (i < m1.numRows) {
        var j = 0
        while (j < m1.numCols) {
          s.maddTo(m1.data, i + m1.numRows*j, m2.data, j + m2.numRows*i, ret, 0)
          j += 1
        }
        i += 1
      }
    }
  }
}



// --------------------------------------
// Lapack operations


trait DenseLapackImplicits {
  implicit def toDenseRealLapackOps   [S <: Scalar.RealTyp]   (m: Dense[S]) = new DenseRealLapackOps(m)
  implicit def toDenseComplexLapackOps[S <: Scalar.ComplexTyp](m: Dense[S]) = new DenseComplexLapackOps(m)
}


object DenseLapackOps {
    /** A \ V -> X where A X ~= V */
  def QRSolve[S <: Scalar](X: Dense[S], A: Dense[S], V: Dense[S], transpose: Boolean)(implicit mb: MatrixBuilder[S, Dense]) {
    require(
        A.numCols == X.numRows && A.numRows == V.numRows && X.numCols == V.numCols, 
        "Cannot divide matrices: [%d, %d] \\ [%d, %d] -> [%d, %d].".format(
            A.numRows, A.numCols, V.numRows, V.numCols, X.numRows, X.numCols 
        )
    )
    require(X.netlib != null, "Netlib library required for division operation.")
    
    val nrhs = V.numCols;
    
    // allocate temporary solution matrix
    // TODO: is X big enough to use directly?
    val Xtmp = mb.zeros(math.max(A.numRows, A.numCols), nrhs)
    val M = if (!transpose) A.numRows else A.numCols;
    for (j <- 0 until nrhs; i <- 0 until M) { Xtmp(i, j) = V(i, j) }

    val newData = A.duplicate(mb);

    // query optimal workspace
    val queryWork = mb.zeros(1, 1)
    val queryInfo = new com.sun.jna.ptr.IntByReference(0);
    X.netlib.gels(
      if (!transpose) "N" else "T",
      A.numRows, A.numCols, nrhs,
      newData.data.buffer, math.max(1,A.numRows),
      Xtmp.data.buffer, math.max(1,math.max(A.numRows,A.numCols)),
      queryWork.data.buffer, -1, queryInfo)
    
    // allocate workspace
    val workSize =
      if (queryInfo.getValue != 0)
        math.max(1, math.min(A.numRows, A.numCols) + math.max(math.min(A.numRows, A.numCols), nrhs));
      else
        math.max(1, X.netlib.readBufferInt(queryWork.data.buffer, 0));
    val work = mb.zeros(workSize, 1)

    // compute factorization
    X.netlib.gels(
      if (!transpose) "N" else "T",
      A.numRows, A.numCols, nrhs,
      newData.data.buffer, math.max(1,A.numRows),
      Xtmp.data.buffer, math.max(1,math.max(A.numRows,A.numCols)),
      work.data.buffer, workSize, queryInfo);

    if (queryInfo.getValue< 0)
      throw new IllegalArgumentException;

    // extract solution
    val N = if (!transpose) A.numCols else A.numRows;
    for (j <- 0 until nrhs; i <- 0 until N) X(i, j) = Xtmp(i, j)
  }
}


class DenseLapackOps[S <: Scalar](self: Dense[S]) {
  def \ (V: Dense[S])(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(V.numRows, V.numCols)
    DenseLapackOps.QRSolve(ret, self, V, false)
    ret
  }
}

class DenseRealLapackOps[S <: Scalar.RealTyp](self: Dense[S]) extends DenseLapackOps(self) {
  def eig(implicit mb: MatrixBuilder[S, Dense]): (Dense[S], Dense[S], Dense[S]) = {
    self.netlib match {
      case netlib: NetlibReal[_] => {
        require(self.numRows == self.numCols, "Cannot find eigenvectors of non-square matrix [%d, %d].".format(self.numRows, self.numCols))
        
        val n = self.numRows

        // Allocate space for the decomposition
        val Wr = mb.zeros(n, 1)
        val Wi = mb.zeros(n, 1)
        val Vr = mb.zeros(n, n)
        
        // Find the needed workspace
        val worksize = mb.zeros(1, 1)
        val info = new com.sun.jna.ptr.IntByReference(0)
        netlib.geev(
          "N", "V", // compute right eigenvectors only
          n,
          netlib.emptyBuffer, math.max(1,n),
          netlib.emptyBuffer, netlib.emptyBuffer,
          netlib.emptyBuffer, math.max(1,n),
          netlib.emptyBuffer, math.max(1,n),
          worksize.data.buffer, -1, info)

        // Allocate the workspace
        val lwork: Int =
          if (info.getValue != 0)
            math.max(1, 4*n)
          else
            math.max(1, netlib.readBufferInt(worksize.data.buffer, 0))
        val work = mb.zeros(lwork, 1)

        // Factor it!
        netlib.geev(
          "N", "V", n,
          self.duplicate.data.buffer, math.max(1,n),
          Wr.data.buffer, Wi.data.buffer,
          netlib.emptyBuffer, math.max(1, n),
          Vr.data.buffer, math.max(1,n),
          work.data.buffer, lwork, info)

        require(info.getValue >= 0, "Error in dgeev argument %d.".format(-info.getValue))
        require(info.getValue <= 0, "Not converged dgeev; only %d of %d eigenvalues computed.".format(info.getValue, self.numRows))
        
        (Wr, Wi, Vr)
      }
      case _ => sys.error("Netlib sgeev/dgeev operation unavailable.")
    }
  }
}

class DenseComplexLapackOps[S <: Scalar.ComplexTyp](self: Dense[S]) extends DenseLapackOps(self) {
  def eig(implicit mb: MatrixBuilder[S, Dense]): (Dense[S], Dense[S]) = {
    self.netlib match {
      case netlib: NetlibComplex[_] => {
        require(self.numRows == self.numCols, "Cannot find eigenvectors of non-square matrix [%d, %d].".format(self.numRows, self.numCols))
        
        val n = self.numRows

        // Allocate space for the decomposition
        val W = mb.zeros(n, 1)
        val Vr = mb.zeros(n, n)

        // Find the needed workspace
        val worksize = mb.zeros(1, 1)
        val info = new com.sun.jna.ptr.IntByReference(0)
        netlib.geev(
          "N", "V", // compute right eigenvectors only
          n,
          netlib.emptyBuffer, math.max(1,n),
          netlib.emptyBuffer,
          netlib.emptyBuffer, math.max(1,n),
          netlib.emptyBuffer, math.max(1,n),
          worksize.data.buffer, -1, netlib.emptyBuffer, info)

        // Allocate the workspace
        val lwork: Int =
          if (info.getValue != 0)
            math.max(1, 4*n)
          else
            math.max(1, netlib.readBufferInt(worksize.data.buffer, 0))
        val work = mb.zeros(lwork, 1)
        val rwork = mb.zeros(n, 1) // 2*n float/double values

        // Factor it!
        netlib.geev(
          "N", "V", n,
          self.duplicate.data.buffer, math.max(1,n),
          W.data.buffer,
          netlib.emptyBuffer, math.max(1, n),
          Vr.data.buffer, math.max(1,n),
          work.data.buffer, lwork, rwork.data.buffer, info)

        require(info.getValue >= 0, "Error in dgeev argument %d.".format(-info.getValue))
        require(info.getValue <= 0, "Not converged dgeev; only %d of %d eigenvalues computed.".format(info.getValue, self.numRows))
        
        (W, Vr)
      }
      case _ => sys.error("Netlib cgeev/zgeev operation unavailable.")
    }
  }
}
