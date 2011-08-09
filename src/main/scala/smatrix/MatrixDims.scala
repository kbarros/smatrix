package smatrix


object MatrixDims {
  def checkDims(numRows: Int, numCols: Int) {
    require(numRows > 0 && numCols > 0,
        "Cannot build matrix with non-positive dimensions [%d, %d].".format(numRows, numCols))
  }
  
  def checkKey(m: MatrixDims, i: Int, j: Int) {
    require(0 <= i && i < m.numRows && 0 <= j && j < m.numCols,
        "Matrix indices out of bounds: [%d %d](%d %d).".format(m.numRows, m.numCols, i, j))
  }
  
  def checkAdd(m1: MatrixDims, m2: MatrixDims) {
    require(m1.numRows == m2.numRows && m1.numCols == m2.numCols,
        "Cannot add/subtract matrices of shape [%d, %d] +- [%d, %d].".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }
  def checkAddTo(m: MatrixDims, ret: MatrixDims) {
    require(m.numRows == ret.numRows && m.numCols == ret.numCols,
        "Cannot add/subtract matrices of shape +-[%d, %d] + [%d, %d] -> [%d, %d].".format(m.numRows, m.numCols, ret.numRows, ret.numCols, ret.numRows, ret.numCols))
  }

  def checkMul(m1: MatrixDims, m2: MatrixDims) {
    require(m1.numCols == m2.numRows,
        "Cannot multiply matrices of shape [%d, %d] * [%d, %d].".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }
  def checkMulTo(m1: MatrixDims, m2: MatrixDims, ret: MatrixDims) {
    checkMul(m1, m2)
    require(
        ret.numRows == m1.numRows &&
        m1.numCols == m2.numRows &&
        m2.numCols == ret.numCols, "Cannot multiply matrices of shape: [%d, %d] * [%d, %d] -> [%d, %d].".format(
            m1.numRows, m1.numCols, m2.numRows, m2.numCols, ret.numRows, ret.numCols))
  }
  
  def checkDot(m1: MatrixDims, m2: MatrixDims) {
    checkMul(m1, m2)
    require(m1.numRows == 1 && m2.numCols == 1,
        "Dot product expects row and column vectors, found [%d, %d] * [%d, %d].".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }
}

trait MatrixDims {
  def numRows: Int
  def numCols: Int
}
