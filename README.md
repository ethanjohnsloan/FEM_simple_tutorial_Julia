# FEM_simple_tutorial_Julia

Performs extremely simple examples of FEM for tutorial purposes. This is not meant to be used for any research applications, only as a beginner learning resource. Supported by Julia 1.12.6. 

This resource accompanies The Finite Element Method: Theory, Implementation, and Applications by M.G. Larson ad F. Berngzon (Springer, 2013).

The "1D" folder has numerous ways of solving example 2.4.2 in the book. The first is, roughly, a translation of Larson and Berngzon's MatLab code. It uses LinearSolve.jl. Another example is the same, except it implements PETSc to solve the equation Ax=b. This is not entirely useful on its own, but familiarity with PETSc is importent in scientific computing and is more condicive to MPI implementation. The MPI implementation here is incomplete.
