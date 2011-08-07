smatrix: *fast matrices in scala*
=================================

A Scala library for high performance matrix computation, with support for native BLAS and LAPACK libraries.


Features
--------

* Optimized for large dense and sparse matrices (no special support for small matrices nor general tensors).

* Fully integrated complex support. Matrix elements, real or complex, are packed in a JVM array or a native buffer (planned). Heap allocation of matrix elements is avoided when possible.

* The only external dependency is [JNA](https://github.com/twall/jna), which is used to call native BLAS and LAPACK libraries. Currently only `vecLib` on OS X is supported.

* Typeclass based design enables user extension, including new matrix representations, matrix operations, and scalar types.

* Natural syntax similar to that of Matlab. In-place matrix operations are also available for performance critical code (planned).


Installation
------------

Distribution is currently source only. Get the code from Github. If you want to make smatrix a dependency for your SBT 0.10 project, modify your `Project` definition in the `build/` directory as follows

    object MyBuild extends Build {
      ...
      lazy val smatrixProject = RootProject( uri("git://github.com/kbarros/smatrix.git") )
      // alternatively, RootProject( file("...") )
      
      lazy val myProject = Project( ... ) dependsOn smatrixProject
    }

For more details on configuration, refer to the [SBT Wiki](https://github.com/harrah/xsbt/wiki/Full-Configuration).

Currently smatrix supports the OS X `vecLib` library for BLAS and LAPACK operations. Support for the open source [ATLAS library](http://math-atlas.sourceforge.net/) is planned. Contact me if you have a specific interest in an alternate library (e.g., Intel MKL).


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
    
    val n = 1000000
    val m = sparse(n, n)                       // Create a 1Mx1M sparse matrix backed by a HashMap
    m(1, 2) = 3                                // Fill in some elements
    m(3, 4) = 6 
    val x = tabulate(n, 1)((i, j) => i)        // Create a dense column vector (x_ij = i)
    val y = m * x                              // Sparse-Dense matrix multiplication


A Work in Progress
------------------

I'm developing smatrix for use in my physics research. The interfaces are subject to change. Your suggestions and contributions are welcome.

**Notes to self, TODO list**

These are things I'm loosely planning to do, but I probably won't do them until there's a need.

- Packed sparse for Sparse*Dense performance
- Low priority implicits? (e.g., general Matrix + Matrix -> Dense operation)
- Expose in place operations
- Get specialization working with parallel "bare" heirarchy
- Add ATLAS backend; JNA name remapping?
- Test on Linux
- Hermitian matrices, perhaps others (Triangular, Tridiagonal, Unitary)
- Matrix solvers for sparse matrices (GMRES, BiCGSTAB)
- Native buffers (GPU implementation?)
