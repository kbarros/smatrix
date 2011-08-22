package smatrix


object Constructors {
  object realFlt extends Constructors[Scalar.RealFlt]
  object realDbl extends Constructors[Scalar.RealDbl]
  object complexFlt extends Constructors[Scalar.ComplexFlt] {
    implicit def floatToComplexf[T <% Float](re: T) = Complexf(re, 0)
    val I = Complexf(0, 1)
  }
  object complexDbl extends Constructors[Scalar.ComplexDbl] {
    implicit def doubleToComplexd[T <% Double](re: T) = Complexd(re, 0)
    val I = Complexd(0, 1)
  }
}


class Constructors[S <: Scalar] {
  // --------------------------------------
  // Dense constructors
  
  def dense(numRows: Int, numCols: Int)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    mb.zeros(numRows, numCols)
  }
  
  def col(elems: S#A*)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val m = dense(elems.size, 1)
    for (i <- 0 until m.numRows) m(i, 0) = elems(i)
    m
  }

  def row(elems: S#A*)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val m = dense(1, elems.size)
    for (j <- 0 until m.numCols) m(0, j) = elems(j)
    m
  }

  def fromRows(row1: Dense[S], rows: Dense[S]*)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    require(row1.numRows == 1 && rows.forall(_.numRows == 1))
    require(rows.forall(_.numCols == row1.numCols))
    val ret = dense(1 + rows.size, row1.numCols)
    ret(0, ::) = row1
    for (i <- rows.indices)
      ret(i+1, ::) = rows(i)
      ret
  }
  
  // --------------------------------------
  // Sparse constructors
  
  def sparse(numRows: Int, numCols: Int)(implicit mb: MatrixBuilder[S, HashSparse]): HashSparse[S] = {
    mb.zeros(numRows, numCols)
  }
  
  def eye(numRows: Int)(implicit mb: MatrixBuilder[S, HashSparse]): HashSparse[S] = {
    val m = sparse(numRows, numRows)
    for (i <- 0 until numRows) m(i, i) = m.scalar.one
    m
  }
}

