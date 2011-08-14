package smatrix


object MatrixDims {
  def checkDims(numRows: Int, numCols: Int) {
    require(numRows > 0 && numCols > 0,
        "Matrix dimensions must be positive, found [%d, %d].".format(numRows, numCols))
  }
  
  def checkKey(m: MatrixDims, i: Int, j: Int) {
    require(0 <= i && i < m.numRows && 0 <= j && j < m.numCols,
        "Matrix indices out of bounds: [%d %d](%d %d).".format(m.numRows, m.numCols, i, j))
  }
  
  def checkAssign(m1: MatrixDims, m2: MatrixDims) {
    require(m1.numRows == m2.numRows && m1.numCols == m2.numCols,
        "Cannot assign matrices of shape [%d, %d] -> [%d, %d].".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }

  def checkAdd(m1: MatrixDims, m2: MatrixDims) {
    require(m1.numRows == m2.numRows && m1.numCols == m2.numCols,
        "Cannot add/subtract matrices of shape [%d, %d] +- [%d, %d].".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }

  def checkMul(m1: MatrixDims, m2: MatrixDims) {
    require(m1.numCols == m2.numRows,
        "Cannot multiply matrices of shape [%d, %d] * [%d, %d].".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }
  
  def checkMulTo(m1: MatrixDims, m2: MatrixDims, ret: MatrixDims) {
    require(
        ret.numRows == m1.numRows &&
        m1.numCols == m2.numRows &&
        m2.numCols == ret.numCols, "Cannot multiply matrices of shape: [%d, %d] * [%d, %d] -> [%d, %d].".format(
            m1.numRows, m1.numCols, m2.numRows, m2.numCols, ret.numRows, ret.numCols))
    require((m1 ne ret) && (m2 ne ret), "Illegal aliasing in matrix product.")
  }
  
  def checkDot(m1: MatrixDims, m2: MatrixDims) {
    checkMul(m1, m2)
    require(m1.numRows == 1 && m2.numCols == 1,
        "Dot product expects row and column vectors, found [%d, %d] * [%d, %d].".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }
  
  def checkTrace(m: MatrixDims) {
    require(m.numRows == m.numCols,
        "Matrix trace expects a square matrix, found [%d, %d].".format(m.numRows, m.numCols))
  }

}

trait MatrixDims {
  def numRows: Int
  def numCols: Int
}
