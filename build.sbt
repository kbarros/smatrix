name := "smatrix"

version := "0.1"

scalaVersion := "2.9.1.RC1"

libraryDependencies ++= Seq(
  "net.java.dev.jna" % "jna" % "3.3.0"
)

resolvers ++= Seq(
  "download.java.net" at "http://download.java.net/maven/2"
)


scalacOptions ++= Seq(
  "-deprecation", "-unchecked"
)
