package smatrix


object Test extends App {
  import Constructors.complexDbl._
  val m = dense(4, 4)
  m(0, 0) = 1

  println((m+m).map(_.re))
  
  val n = 500
  val m3 = tabulate(n, n) { case (i, j) => i + 2*j }
  val m4 = m3 + m3.tran
  val x = tabulate(n, 1) { case (i, j) => i }
  val (v, w) = m3.eig
  println(m3 * w(::,0) / v(0) - w(::,0))
}
