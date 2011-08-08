package smatrix


object Test extends App {
  import Constructors.complexDbl._
  val m = dense(4, 4)

  // types should be inferred here
  val m2 = (m+m).map(_.re).map(_ + I)
  
  val n = 500
  val m3 = tabulate(n, n) { case (i, j) => i + 2*j }
  val m4 = m3 + m3.tran
  val x = tabulate(n, 1) { case (i, j) => i }
  val (v, w) = m3.eig
  println(m3 * w(::,0) / v(0) - w(::,0))
  println("should be zero\n\n")
  
  val s = sparse(2, 4)
  s(0, 0) = 2+I
  s(1, 3) = 2-I
  println(s * col(1, 1, 1, 2) - col(2+I, 4-2*I))
  println("should be zero\n\n")
  
  val s2 = sparse(2, 2)
  s + s2
}
