package smatrix

import java.nio.{DoubleBuffer, FloatBuffer}
import java.nio.Buffer


object Scalar {
  trait RealTyp    extends Scalar
  trait ComplexTyp extends Scalar
  
  trait RealDbl    extends RealTyp    { type A = Double;   type Raw = Double; type Buf = DoubleBuffer }
  trait RealFlt    extends RealTyp    { type A = Float;    type Raw = Float;  type Buf = FloatBuffer }
  trait ComplexDbl extends ComplexTyp { type A = Complexd; type Raw = Double; type Buf = DoubleBuffer }
  trait ComplexFlt extends ComplexTyp { type A = Complexf; type Raw = Float;  type Buf = FloatBuffer }
}

trait Scalar {
  type A
  type Raw
  type Buf <: Buffer
}


object ScalarOps {
  implicit object RealFlt extends RealFlt
  implicit object RealDbl extends RealDbl
  implicit object ComplexFlt extends ComplexFlt
  implicit object ComplexDbl extends ComplexDbl

  trait RealFlt extends ScalarOps[Scalar.RealFlt] {
    def add(a: Float, b: Float): Float = a + b
    def sub(a: Float, b: Float): Float = a - b
    def mul(a: Float, b: Float): Float = a * b
    def div(a: Float, b: Float): Float = a / b
    def neg(a: Float): Float = -a
    def conj(a: Float): Float = a
    def zero: Float = 0.0f
    def one: Float = 1.0f
    
    def components = 1
    def read(data: Data, i: Int): Float = data(i)
    def write(data: Data, i: Int, x: Float): Unit = data(i) = x
    override def maddTo(d1: Data, i1: Int, x2: Float, d3: Data, i3: Int) {
      d3(i3) += d1(i1) * x2
    }
    override def maddTo(d1: Data, i1: Int, d2: Data, i2: Int, d3: Data, i3: Int) {
      d3(i3) += d1(i1) * d2(i2)
    }
  }

  trait RealDbl extends ScalarOps[Scalar.RealDbl] {
    def add(a: Double, b: Double): Double = a + b
    def sub(a: Double, b: Double): Double = a - b
    def mul(a: Double, b: Double): Double = a * b
    def div(a: Double, b: Double): Double = a / b
    def neg(a: Double): Double = -a
    def conj(a: Double): Double = a
    def zero: Double = 0.0
    def one: Double = 1.0
    
    def components = 1
    def read(data: Data, i: Int): Double = data(i)
    def write(data: Data, i: Int, x: Double): Unit = data(i) = x
    override def maddTo(d1: Data, i1: Int, x2: Double, d3: Data, i3: Int) {
      d3(i3) += d1(i1) * x2
    }
    override def maddTo(d1: Data, i1: Int, d2: Data, i2: Int, d3: Data, i3: Int) {
      d3(i3) += d1(i1) * d2(i2)
    }
  }

  trait ComplexFlt extends ScalarOps[Scalar.ComplexFlt] {
    def add(a: Complexf, b: Complexf): Complexf = a + b
    def sub(a: Complexf, b: Complexf): Complexf = a - b
    def mul(a: Complexf, b: Complexf): Complexf = a * b
    def div(a: Complexf, b: Complexf): Complexf = a / b
    def neg(a: Complexf): Complexf = -a
    def conj(a: Complexf): Complexf = a.conj
    def zero: Complexf = Complexf(0f, 0f)
    def one: Complexf = Complexf(1f, 0f)
    
    def components = 2
    def read(data: Data, i: Int) = Complexf(data(2*i+0), data(2*i+1))
    def write(data: Data, i: Int, x: Complexf) {
      data(2*i+0) = x.re
      data(2*i+1) = x.im
    }
    override def maddTo(d1: Data, i1: Int, x2: Complexf, d3: Data, i3: Int) {
      val x1_re = d1(2*i1+0)
      val x1_im = d1(2*i1+1)
      d3(2*i3+0) += x1_re*x2.re - x1_im*x2.im
      d3(2*i3+1) += x1_re*x2.im + x1_im*x2.re
    }
    override def maddTo(d1: Data, i1: Int, d2: Data, i2: Int, d3: Data, i3: Int) {
      val x1_re = d1(2*i1+0)
      val x1_im = d1(2*i1+1)
      val x2_re = d2(2*i2+0)
      val x2_im = d2(2*i2+1)
      d3(2*i3+0) += x1_re*x2_re - x1_im*x2_im
      d3(2*i3+1) += x1_re*x2_im + x1_im*x2_re
    }
  }

  trait ComplexDbl extends ScalarOps[Scalar.ComplexDbl] {
    def add(a: Complexd, b: Complexd): Complexd = a + b
    def sub(a: Complexd, b: Complexd): Complexd = a - b
    def mul(a: Complexd, b: Complexd): Complexd = a * b
    def div(a: Complexd, b: Complexd): Complexd = a / b
    def neg(a: Complexd): Complexd = -a
    def conj(a: Complexd): Complexd = a.conj
    def zero: Complexd = Complexd(0d, 0d)
    def one: Complexd = Complexd(1d, 0d)
    
    def components = 2
    def read(data: Data, i: Int) = Complexd(data(2*i+0), data(2*i+1))
    def write(data: Data, i: Int, x: Complexd) {
      data(2*i+0) = x.re
      data(2*i+1) = x.im
    }
    override def maddTo(d1: Data, i1: Int, x2: Complexd, d3: Data, i3: Int) {
      val x1_re = d1(2*i1+0)
      val x1_im = d1(2*i1+1)
      d3(2*i3+0) += x1_re*x2.re - x1_im*x2.im
      d3(2*i3+1) += x1_re*x2.im + x1_im*x2.re
    }
    override def maddTo(d1: Data, i1: Int, d2: Data, i2: Int, d3: Data, i3: Int) {
      val x1_re = d1(2*i1+0)
      val x1_im = d1(2*i1+1)
      val x2_re = d2(2*i2+0)
      val x2_im = d2(2*i2+1)
      d3(2*i3+0) += x1_re*x2_re - x1_im*x2_im
      d3(2*i3+1) += x1_re*x2_im + x1_im*x2_re
    }
  }
}


trait GenScalarOps[@specialized(Float, Double) A, @specialized(Float, Double) Raw, Buf <: Buffer] {
  def add(a: A, b: A): A
  def sub(a: A, b: A): A
  def mul(a: A, b: A): A
  def div(a: A, b: A): A
  def neg(a: A): A
  def conj(a: A): A
  def zero: A
  def one: A
  
  type Data = RawData[Raw, Buf]
  def components: Int
  def read(data: Data, i: Int): A
  def write(data: Data, i: Int, x: A)

  def maddTo(d1: Data, i1: Int, x2: A, d3: Data, i3: Int) {
    val x1 = read(d1, i1)
    val x3 = read(d3, i3)
    write(d3, i3, add(mul(x1, x2), x3))
  }
  
  def maddTo(d1: Data, i1: Int, d2: Data, i2: Int, d3: Data, i3: Int) {
    val x1 = read(d1, i1)
    val x2 = read(d2, i2)
    val x3 = read(d3, i3)
    write(d3, i3, add(mul(x1, x2), x3))
  }
}
trait ScalarOps[S <: Scalar] extends GenScalarOps[S#A, S#Raw, S#Buf]
