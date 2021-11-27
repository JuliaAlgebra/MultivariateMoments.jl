# Atoms extration

## Vectorized matrix

```@docs
SymMatrix
VectorizedHermitianMatrix
square_getindex
symmetric_setindex!
```

## Moment matrix

```@docs
MomentMatrix
moment_matrix
```

## Atomic measure

```@docs
WeightedDiracMeasure
AtomicMeasure
```

## Atoms extraction

Given a `MomentMatrix` with a positive semidefinite moment matrix,
an algebraic system for which the set of solution is a superset of the support of the measure.
If the measure is atomic and the `MomentMatrix` contains enough moments,
the algebraic system will only have a finite number of solutions which are the centers
of the diracs of the measure.

```@docs
MultivariateMoments.computesupport!
extractatoms
```

This system is obtained from a low rank cholesky decomposition of the moment matrix.
This decomposition can either be obtained by a cholesky or SVD decomposition from which we remove the rows corresponding to the negligeable eigenvalues/singular values.

```@docs
LowRankChol
ShiftChol
SVDChol
MultivariateMoments.lowrankchol
```

Once the center of the atoms are determined, a linear system is solved to determine
the weights corresponding to each dirac.
By default, [`MomentMatrixWeightSolver`](@ref) is used by [`extractatoms`](@ref) so that if there are small differences between moment values corresponding to the same monomial in the matrix
(which can happen if these moments were computed numerically by a semidefinite proramming solvers, e.g., with [SumOfSquares](https://github.com/jump-dev/SumOfSquares.jl)),
the linear system handles that automatically.
```@docs
MomentMatrixWeightSolver
MomentVectorWeightSolver
```
