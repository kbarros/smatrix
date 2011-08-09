package smatrix


import java.nio.FloatBuffer
import java.nio.DoubleBuffer
import java.nio.Buffer


object RawData {
  class FltArray(val size: Int) extends RawData[Float, FloatBuffer] {
    val array = new Array[Float](size)
    val buffer = java.nio.FloatBuffer.wrap(array)
    def apply(i: Int): Float = array(i)
    def update(i: Int, x: Float) { array(i) = x }
  }

  class DblArray(val size: Int) extends RawData[Double, DoubleBuffer] {
    val array = new Array[Double](size)
    val buffer = java.nio.DoubleBuffer.wrap(array)
    def apply(i: Int): Double = array(i)
    def update(i: Int, x: Double): Unit = array(i) = x
  }
}


trait RawData[@specialized(Float, Double) Raw, Buf <: Buffer] {
  def size: Int
  def buffer: Buf
  def apply(i: Int): Raw
  def update(i: Int, x: Raw)
  def copyTo(that: RawData[Raw, Buf]) { for (i <- 0 until size) that(i) = this(i) }
  def dispose() {}
}

