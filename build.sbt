name := "smatrix"

version := "0.1"

scalaVersion := "2.10.3"

libraryDependencies ++= Seq(
  "net.java.dev.jna" % "jna" % "3.3.0"
)

resolvers ++= Seq(
  "download.java.net" at "http://download.java.net/maven/2"
)


scalacOptions ++= Seq(
  "-deprecation", "-unchecked",
  "-language:implicitConversions", "-language:higherKinds"
)
