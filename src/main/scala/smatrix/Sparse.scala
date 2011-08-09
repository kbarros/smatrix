package smatrix

import collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer


object Sparse extends SparseBuilders with SparseAdders with SparseMultipliers


abstract class Sparse[S <: Scalar : ScalarOps, +Repr[s <: Scalar] <: Sparse[s, Repr]]
    (numRows: Int, numCols: Int) extends Matrix[S, Repr](numRows, numCols) { self: Repr[S] =>

  override def transform(f: S#A => S#A): this.type = {
    require(f(scalar.zero) == scalar.zero, "Sparse transformation function must preserve zero.")
    for ((i, j) <- definedIndices) { this(i, j) = f(this(i, j)) }
    this
  }

  override def map[A2, S2 <: Scalar{type A=A2}, That[s <: Scalar] >: Repr[s] <: Matrix[s, That]]
      (f: S#A => S2#A)(implicit scalar2: ScalarOps[S2], mb: MatrixBuilder[S2, That]): That[S2] = {
    require(f(scalar.zero) == scalar2.zero, "Sparse map function must preserve zero.")
    mb.map(this)(f)
  }
}


class HashSparse[S <: Scalar : ScalarOps](numRows: Int, numCols: Int)
    extends Sparse[S, HashSparse](numRows, numCols) {
  val data = new HashMap[(Int, Int), S#A]()
  override val description = "Sparse Hashtable"
  override def definedIndices: Iterable[(Int, Int)] = data.keys
  override def apply(i: Int, j: Int): S#A = data.getOrElse((i, j), scalar.zero)
  override def update(i: Int, j: Int, x: S#A) { data((i, j)) = x }
  
  def toPacked(implicit db: ScalarBuilder[S]): PackedSparse[S] = {
    val ret = PackedSparse.fromIndices(numRows, numCols, definedIndices)
    for ((i, j) <- definedIndices) {
      ret(i, j) = this(i, j)
    }
    ret
  }
}


object PackedSparse {
  def fromIndices[S <: Scalar : ScalarOps : ScalarBuilder]
      (numRows: Int, numCols: Int, indices: Iterable[(Int, Int)]) : PackedSparse[S] = {
    new PackedSparse[S](numRows, numCols) {
      override val data = implicitly[ScalarBuilder[S]].build(indices.size)
      override val definedCols = {
        val cols = Array.fill(numRows)(new ArrayBuffer[Int]())
        for ((i, j) <- indices) {
          cols(i).append(j)
        }
        cols.map(_.sorted.toArray)
      }
      override val definedColsAccum = {
        val ret = Array.fill(numRows)(0)
        var acc = 0
        for (i <- 0 until numRows) {
          ret(i) = acc
          acc += definedCols(i).size
        }
        ret
      }
    }
  }
}

abstract class PackedSparse[S <: Scalar : ScalarOps]
    (numRows: Int, numCols: Int) extends Sparse[S, PackedSparse](numRows, numCols) {
  val data: RawData[S#Raw, S#Buf]
  val definedCols: Array[Array[Int]]
  val definedColsAccum: Array[Int]
  
  override val description = "Sparse Packed"
  override def definedIndices: Iterable[(Int, Int)] = {
    for (i <- 0 until numRows; j <- definedCols(i)) yield (i, j)
  }
  override def apply(i: Int, j: Int): S#A = {
    indexOption(i, j) match {
      case None => scalar.zero
      case Some(idx) => scalar.read(data, idx)
    }
  }
  override def update(i: Int, j: Int, x: S#A) {
    indexOption(i, j) match {
      case None => require(false, "Cannot write to undefined index (%d %d) in PackedSparse matrix.".format(i, j))
      case Some(idx) => scalar.write(data, idx, x)
    }
  }
  
  def indexOption(i: Int, j: Int): Option[Int] = {
    MatrixDims.checkKey(this, i, j)
    val idx = definedCols(i).indexOf(j)
    if (idx == -1) None else Some(idx + definedColsAccum(i))
  }
}



trait SparseBuilders {
  implicit def hashSparseBuilder[S <: Scalar : ScalarOps] = new MatrixBuilder[S, HashSparse] {
    def zeros(numRows: Int, numCols: Int): HashSparse[S] = {
      new HashSparse[S](numRows, numCols)
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

  implicit def packedSparseBuilder[S <: Scalar : ScalarOps : ScalarBuilder] = new MatrixBuilder[S, PackedSparse] {
    def zeros(numRows: Int, numCols: Int): PackedSparse[S] = {
      PackedSparse.fromIndices(numRows, numCols, Seq())
    }
    def duplicate(m: PackedSparse[S]): PackedSparse[S] = {
      val ret = PackedSparse.fromIndices(m.numRows, m.numCols, m.definedIndices)
      for ((i, j) <- m.definedIndices) { ret(i, j) = m(i, j) }
      ret
    }
    def transpose(m: PackedSparse[S]): PackedSparse[S] = {
      val indices = m.definedIndices.map { case (i, j) => (j, i) }
      val ret = PackedSparse.fromIndices(m.numCols, m.numRows, indices)
      for ((i, j) <- indices) { ret(i, j) = m(j, i) }
      ret
    }
    def map[S0 <: Scalar](m: PackedSparse[S0])(f: S0#A => S#A): PackedSparse[S] = {
      val ret = PackedSparse.fromIndices(m.numRows, m.numCols, m.definedIndices)
      for ((i, j) <- m.definedIndices) { ret(i, j) = f(m(i, j)) }
      ret
    }
  }
}


trait SparseAdders {
  implicit def sparseIntoSparseAdder[S <: Scalar, M[s <: Scalar] <: Sparse[s, M]] = new MatrixSingleAdder[S, M, HashSparse] {
    def addTo(neg: Boolean, m: M[S], ret: HashSparse[S]) = {
      MatrixDims.checkAddTo(m, ret)
      for ((i, j) <- m.definedIndices)
        ret(i, j) = if (neg) ret.scalar.sub(ret(i, j), m(i, j)) else ret.scalar.add(ret(i, j), m(i, j))
    }
  }
  
  implicit def sparseIntoDenseAdder[S <: Scalar, M[s <: Scalar] <: Sparse[s, M]] = new MatrixSingleAdder[S, M, Dense] {
    def addTo(neg: Boolean, m: M[S], ret: Dense[S]) = {
      MatrixDims.checkAddTo(m, ret)
      for ((i, j) <- m.definedIndices)
        ret(i, j) = if (neg) ret.scalar.sub(ret(i, j), m(i, j)) else ret.scalar.add(ret(i, j), m(i, j))
    }
  }

  implicit def sparseSparseAdder[S <: Scalar, M1[s <: Scalar] <: Sparse[s, M1], M2[s <: Scalar] <: Sparse[s, M2]] =
    new MatrixAdder[S, M1, M2, HashSparse](sparseIntoSparseAdder, sparseIntoSparseAdder)

  implicit def denseSparseAdder[S <: Scalar, M2[s <: Scalar] <: Sparse[s, M2]] =
    new MatrixAdder[S, Dense, M2, Dense](Dense.denseIntoDenseAdder, sparseIntoDenseAdder)

  implicit def sparseDenseAdder[S <: Scalar, M1[s <: Scalar] <: Sparse[s, M1]] =
    new MatrixAdder[S, M1, Dense, Dense](sparseIntoDenseAdder, Dense.denseIntoDenseAdder)
}


trait SparseMultipliersLowPriority {
  implicit def sparseDenseMultiplier[S <: Scalar, M[s <: Scalar] <: Sparse[s, M]] = new MatrixMultiplier[S, M, Dense, Dense] {
    def maddTo(m1: M[S], m2: Dense[S], ret: Dense[S]) {
      MatrixDims.checkMulTo(m1, m2, ret)
      for ((i, k) <- m1.definedIndices;
           j <- 0 until ret.numCols) {
        ret.scalar.maddTo(m2.data, m2.index(k, j), m1(i, k), ret.data, ret.index(i, j))
      }
    }
  }
  implicit def denseSparseMultiplier[S <: Scalar, M[s <: Scalar] <: Sparse[s, M]] = new MatrixMultiplier[S, Dense, M, Dense] {
    def maddTo(m1: Dense[S], m2: M[S], ret: Dense[S]) {
      MatrixDims.checkMulTo(m1, m2, ret)
      for ((k, j) <- m2.definedIndices;
           i <- 0 until ret.numRows) {
        ret.scalar.maddTo(m1.data, m1.index(i, k), m2(k, j), ret.data, ret.index(i, j))
      }
    }
  }
}
trait SparseMultipliers extends SparseMultipliersLowPriority {
  implicit def packedSparseDenseMultiplier[S <: Scalar] = new MatrixMultiplier[S, PackedSparse, Dense, Dense] {
    def maddTo(m1: PackedSparse[S], m2: Dense[S], ret: Dense[S]) {
      MatrixDims.checkMulTo(m1, m2, ret)
      // ret_ij += \sum_k m1_ik m2_kj
      for (i <- 0 until m1.numRows;
           (k, idx1) <- m1.definedCols(i).zipWithIndex;
           j <- 0 until ret.numCols) {
        ret.scalar.maddTo(m1.data, idx1+m1.definedColsAccum(i), m2.data, m2.index(k, j), ret.data, ret.index(i, j))
      }
    }
  }
}
