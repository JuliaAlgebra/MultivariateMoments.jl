# Atoms extration

## Moment matrix

```@docs
MomentMatrix
moment_matrix
SymMatrix
VectorizedHermitianMatrix
symmetric_setindex!
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
