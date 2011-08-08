package smatrix

import collection.mutable.HashMap


object Sparse extends SparseBuilders with SparseAdders with SparseMultipliers

trait Sparse[S <: Scalar, +Repr[s <: Scalar] <: Sparse[s, Repr]] extends Matrix[S, Repr] { self: Repr[S] =>
  def definedIndices: Iterable[(Int, Int)]

  override def transform(f: S#A => S#A): this.type = {
    for ((i, j) <- definedIndices) { this(i, j) = f(this(i, j)) }
    this
  }

  override def toDense(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(numRows, numCols)
    for ((i, j) <- definedIndices) { ret(i, j) = this(i, j) }
    ret
  }

  override def map[A2, S2 <: Scalar{type A=A2}, That[s <: Scalar] >: Repr[s] <: Matrix[s, That]]
      (f: S#A => S2#A)(implicit scalar2: ScalarOps[S2], mb: MatrixBuilder[S2, That]): That[S2] = {
    require(f(scalar.zero) == scalar2.zero, "Map on sparse matrix must preserve zero")
    mb.map(this)(f)
  }
}


trait HashSparse[S <: Scalar] extends Sparse[S, HashSparse] {
  val data = new HashMap[(Int, Int), S#A]()
  override val description = "Sparse Hashtable"
  override def definedIndices: Iterable[(Int, Int)] = data.keys
  override def apply(i: Int, j: Int): S#A = data.getOrElse((i, j), scalar.zero)
  override def update(i: Int, j: Int, x: S#A) { data((i, j)) = x }
  
  def toPacked(implicit db: RawData.Builder[S#Raw, S#Buf]): PackedSparse[S] = new PackedSparse[S] {
    override val scalar = HashSparse.this.scalar
    override val numRows = HashSparse.this.numRows
    override val numCols = HashSparse.this.numCols
    for (k <- HashSparse.this.data.keys) { data(k) = HashSparse.this.data(k) }
  }
}


trait PackedSparse[S <: Scalar] extends Sparse[S, PackedSparse] {
  val data = new HashMap[(Int, Int), S#A]()
  override val description = "Sparse Packed"
  override def definedIndices: Iterable[(Int, Int)] = data.keys
  override def apply(i: Int, j: Int): S#A = data.getOrElse((i, j), scalar.zero)
  override def update(i: Int, j: Int, x: S#A) { data((i, j)) = x }
}



trait SparseBuilders {
  implicit def hashSparseBuilder[S <: Scalar](implicit so: ScalarOps[S]) = new MatrixBuilder[S, HashSparse] {
    def zeros(numRows: Int, numCols: Int) = {
      MatrixDims.checkDims(numRows, numCols)
      val nr = numRows
      val nc = numCols
      new HashSparse[S] {
        override val scalar = so
        override val numRows = nr
        override val numCols = nc
      }
    }
    
    def duplicate(m: HashSparse[S]): HashSparse[S] = {
      val ret = zeros(m.numRows, m.numCols)
      for (k <- m.data.keys) { ret.data(k) = m.data(k) }
      ret
    }
    
    def transpose(m: HashSparse[S]): HashSparse[S] = {
      val ret = zeros(m.numCols, m.numRows)
      for ((i, j) <- m.data.keys) { ret.data((j, i)) = m.data((i, j)) }
      ret
    }
    
    def map[S0 <: Scalar](m: HashSparse[S0])(f: S0#A => S#A): HashSparse[S] = {
      val ret = zeros(m.numRows, m.numCols)
      for (k <- m.data.keys) { ret.data(k) = f(m.data(k)) }
      ret
    }
  }
}


trait SparseAdders {
  implicit def sparseSparseAdder[S <: Scalar, M[s <: Scalar] <: Sparse[s, M]] = new MatrixAdder[S, M, HashSparse] {
    def addTo(neg: Boolean, m: M[S], ret: HashSparse[S]) = {
      MatrixDims.checkAddTo(m, ret)
      for ((i, j) <- m.definedIndices)
        ret(i, j) = if (neg) ret.scalar.sub(ret(i, j), m(i, j)) else ret.scalar.add(ret(i, j), m(i, j))
    }
  }
}

trait SparseMultipliers {
  implicit def hashSparseDenseMultiplier[S <: Scalar] = new MatrixMultiplier[S, HashSparse, Dense, Dense] {
    def maddTo(m1: HashSparse[S], m2: Dense[S], ret: Dense[S]) {
      MatrixDims.checkMulTo(m1, m2, ret)
      // ret_ij += \sum_k m1_ik m2_kj
      for ((i, k) <- m1.data.keys;
           val m1_ik = m1.apply(i, k); 
           j <- 0 until ret.numCols) {
        ret.scalar.maddTo(m2.data, m2.index(k, j), m1_ik, ret.data, ret.index(i, j))
      }
    }
  }
  implicit def denseHashSparseMultiplier[S <: Scalar] = new MatrixMultiplier[S, Dense, HashSparse, Dense] {
    def maddTo(m1: Dense[S], m2: HashSparse[S], ret: Dense[S]) {
      MatrixDims.checkMulTo(m1, m2, ret)
      // ret_ij += \sum_k m1_ik m2_kj
      for ((k, j) <- m2.data.keys;
           val m2_kj = m1.apply(k, j); 
           i <- 0 until ret.numRows) {
        ret.scalar.maddTo(m1.data, m1.index(i, k), m2_kj, ret.data, ret.index(i, j))
      }
    }
  }
}
