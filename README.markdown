smatrix: *fast matrices in scala*
=================================

A Scala library for high performance matrix computation, with support for native BLAS, LAPACK, and ARPACK libraries.


Features
--------

* Optimized for large dense and sparse matrices (no special support for small matrices nor general tensors).

* Fully integrated complex support. Matrix elements, real or complex, are packed in a JVM array or a native buffer (planned).

* The only external dependency is [JNA](https://github.com/twall/jna), which is used to load netlib libraries.

* Typeclass based design to enable user extension (e.g., matrix representations and operations).

* DSL syntax similar to that of Matlab.


Installation
------------

Source distribution only. To include from SBT, modify your project definition in the `build/` directory:

    object MyBuild extends Build {
      ...
      lazy val smatrixProject = RootProject( uri("git://github.com/kbarros/smatrix.git") )
      // alternatively, RootProject( file("...") )
      
      lazy val myProject = Project( ... ) dependsOn smatrixProject
    }

More details on the [SBT Wiki](https://github.com/harrah/xsbt/wiki/Full-Configuration).

Currently I've hardcoded support for `vecLib`, `acml`, and `arpack` native libraries. Others should be easy to add.


Usage example
--------------

    import smatrix._
    import Constructors.complexDbl._ // other types possible (realDbl, complexFlt, etc.)
    
    // -----------------------
    // Dense Matrices
    
    val m = fromRows(row(1, 2+I), row(3-I, 4)) // 2x2 matrix of type smatrix.Dense[smatrix.Scalar.ComplexDbl]
    val (v, w) = m.eig                         // Get eigenvalues v, eigenvectors w (packed as columns)
    w(::,0) / v(0)                             // First eigenvector, scaled by eigenvalue
    m * w(::,0) / v(0) - w(::,0)               // Zero, up to numerical error
    val x = col(3, 4)                          // A column vector
    m * (m \ x) - x                            // Zero, up to numerical error
    
    // -----------------------
    // Sparse Matrices
    
    val n = 100000
    val m = sparse(n, n)                       // Create a 1Mx1M sparse matrix backed by a HashMap
    m(1, 2) = 3                                // Fill in some elements
    m(3, 4) = 6 
    val x = tabulate(n, 1)((i, j) => i)        // Create a dense column vector (x_ij = i)
    val y = m * x                              // Sparse-Dense matrix multiplication
    val (evals, evecs) =                       // Arpack calculates smallest eigenvalue/vector
      h.toPacked.eig(nev=1, which="SR", tol=1e-4)
    

A Work in Progress
------------------

Still to be done:

- Support for special matrices beyond dense and sparse (e.g., Hermitian, tridiagonal, ...)
- Matrix solvers for sparse matrices (GMRES, BiCGSTAB)
- Native buffers (GPU implementation?)
