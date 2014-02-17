package smatrix


case class Complexd(re: Double, im: Double) {
  def +(that: Double): Complexd = Complexd(re+that, im)
  def -(that: Double): Complexd = Complexd(re-that, im)
  def *(that: Double): Complexd = Complexd(re*that, im*that)
  def /(that: Double): Complexd = Complexd(re/that, im/that)

  def +(that: Complexd): Complexd = Complexd(re+that.re, im+that.im)
  def -(that: Complexd): Complexd = Complexd(re-that.re, im-that.im)
  def *(that: Complexd): Complexd = Complexd(re*that.re - im*that.im, re*that.im + im*that.re)
  def /(that: Complexd): Complexd = this*that.conj / that.abs2

  def exp = {
    val a = math.exp(re)
    Complexd(a * math.cos(im), a * math.sin(im))
  }
  def log = {
    val r2 = re*re + im*im
    val theta = math.atan2(im, re)
    Complexd(0.5*math.log(r2), theta)
  }
  
  def unary_- = Complexd(-re, -im)
  def abs = math.sqrt(abs2)
  def abs2 = re*re + im*im
  def conj = Complexd(re, -im)
  
  def toComplexd = this
  def toComplexf = new Complexf(re.toFloat, im.toFloat)
  
  override def equals(that : Any) = that match {
    case that : Complexd => re == that.re && im == that.im
    case that : Complexf => re == that.re && im == that.im
    case that : Double => re == that && im == 0
    case that : Float => re == that && im == 0
    case that : Long => re == that && im == 0
    case that : Int => re == that && im == 0
    case that : Short => re == that && im == 0
    case _ => false
  }
  
  override def toString = {
    val Tol = 1e-12
    if (math.abs(im) < Tol)
      re.toString
    else if (math.abs(re) < Tol)
      im+"i"
    else {
      if (im >= 0) {
        // re+" + "+im+"i"
        "%g +%gi".format(re, +im)
      }
      else {
        // re+" - "+(-im)+"i"
        "%g -%gi".format(re, -im)
      }
    }
  }
}

case class Complexf(re: Float, im: Float) {
  def +(that: Float): Complexf = Complexf(re+that, im)
  def -(that: Float): Complexf = Complexf(re-that, im)
  def *(that: Float): Complexf = Complexf(re*that, im*that)
  def /(that: Float): Complexf = Complexf(re/that, im/that)

  def +(that: Complexf): Complexf = Complexf(re+that.re, im+that.im)
  def -(that: Complexf): Complexf = Complexf(re-that.re, im-that.im)
  def *(that: Complexf): Complexf = Complexf(re*that.re - im*that.im, re*that.im + im*that.re)
  def /(that: Complexf): Complexf = this*that.conj / that.abs2

  def exp = {
    val a = math.exp(re)
    Complexf((a * math.cos(im)).toFloat, (a * math.sin(im)).toFloat)
  }
  def log = {
    val r2 = re*re + im*im
    val theta = math.atan2(im, re)
    Complexf((0.5*math.log(r2)).toFloat, theta.toFloat)
  }
  
  def unary_- = Complexf(-re, -im)
  def abs = math.sqrt(abs2)
  def abs2 = re*re + im*im
  def conj = Complexf(re, -im)
  
  def toComplexd = new Complexd(re, im)
  def toComplexf = this
  
  override def equals(that : Any) = that match {
    case that : Complexd => re == that.re && im == that.im
    case that : Complexf => re == that.re && im == that.im
    case that : Double => re == that && im == 0
    case that : Float => re == that && im == 0
    case that : Long => re == that && im == 0
    case that : Int => re == that && im == 0
    case that : Short => re == that && im == 0
    case _ => false
  }
  
  override def toString = {
    val Tol = 1e-12
    if (math.abs(im) < Tol)
      re.toString
    else if (math.abs(re) < Tol)
      im+"i"
    else {
      if (im >= 0)
        re+" + "+im+"i"
      else
        re+" - "+(-im)+"i"
    }
  }
}
