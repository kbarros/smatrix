package smatrix

import scala.util.Random


object Test extends App {
  import smatrix._
  import Constructors.complexDbl._
  val rand = new Random(0)

  // basicTests()
  // testArpack()
  testArpackSpeed()
  
  def basicTests() {
    val m = dense(4, 4)

    // types should be inferred here
    val m2 = (m+m).map(_.re).map(_ + I)

    val n = 500
    val m3 = dense(n, n).tabulate { case (i, j) => i + 2*j }
    val m4 = m3 + m3.tran
    val x = dense(n, 1).tabulate { case (i, j) => i }
    val (v, w) = m3.eig
    println(m3 * w(::,0) / v(0) - w(::,0))
    println("should be zero\n\n")

    val s = sparse(2, 4)
    s(0, 0) = 2+I
    s(1, 3) = 2-I
    println(s * col(1, 1, 1, 2) - col(2+I, 4-2*I))
    println("should be zero\n\n")
  }
  
  def randomSparse(n: Int, fillFactor: Int) = {
    val h = sparse(n, n)
    for (iter <- 0 until fillFactor*n) {
      val i = rand.nextInt(n)
      val j = rand.nextInt(n)
      val r = rand.nextGaussian() + rand.nextGaussian() * I
      h(i, j) += r
      h(j, i) += r.conj
    }
    h.toPacked
  }
  
  def testArpack() {
    val h = randomSparse(500, 4)
    
    val (evals0, evecs0) = h.eig(nev=1, which="SR", tol=0)
    val (evals1, evecs1) = h.eig(nev=1, which="LR", tol=0)
    
    val evalsFull = h.toDense.eig._1.toArray.sortBy(_.re)
    
    println(s"Should be zero: ${evalsFull.head - evals0(0)}")
    println(s"Should be zero: ${evalsFull.last - evals1(0)}")
    println(s"Should be zero: ${(h*evecs0 - evecs0*diag(evals0)).norm2}") 
    println(s"Should be zero: ${(h*evecs1 - evecs1*diag(evals1)).norm2}") 
  }
  
  def testArpackSpeed() {
    val n = 50000
    println(s"Building ${n}x${n} matrix")
    val h = randomSparse(n, 4)
    
    print("Calculating min eigenvalue...")
    val t1 = System.currentTimeMillis()
    val (evals0, evecs0) = h.eig(nev=1, which="SR", tol=1e-4)
    val t2 = System.currentTimeMillis()
    println(s"done. Elapsed time ${(t2-t1)/1000.0}s")
    
    println(s"Should be zero: ${(h*evecs0 - evecs0*diag(evals0)).norm2}") 
  }
}
