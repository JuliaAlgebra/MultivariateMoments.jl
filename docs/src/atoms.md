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

Given a `MomentMatrix` containing the moments of an atomic measure,
[`atomic_measure`](@ref) attempts to recover the dirac centers and weights
by first computing an algebraic system with the atom centers as solution
with [`compute_support!`](@ref).

```@docs
compute_support!
atomic_measure
```

For the first step of [`compute_support!`](@ref), there are two approaches.
The first one is to exploit the *flat extension* to directly obtain the multiplication
matrices.

```@docs
FlatExtension
```

The second approach is to first obtain the image space of the moment matrix,
represented as a [`MacaulayNullspace`](@ref)
and to then compute the multiplication matrices from this image space.
This is implemented by the [`ImageSpaceSolver`](@ref).

```@docs
MacaulayNullspace
ImageSpaceSolver
```

This image space is obtained from a low rank LDLT decomposition of the moment matrix.
This decomposition can either be obtained by a cholesky or SVD decomposition from which we remove the rows corresponding to the negligeable eigen/singular values.

```@docs
LowRankLDLTAlgorithm
ShiftCholeskyLDLT
SVDLDLT
low_rank_ldlt
LowRankLDLT
```

The choice of cutoff between the significant and negligeable eigen/singular values is
parametrized by the following interface:
```@docs
RankCheck
rank_from_singular_values
accuracy
doubt
```

The rank check can be chosen among the following:
```@docs
UserRank
FixedRank
FixedRanks
AbsoluteRankTol
LeadingRelativeRankTol
DifferentialRankTol
LargestDifferentialRank
FallbackRank
```

Given the [`MacaulayNullspace`](@ref), there are two approaches implemented
to obtain the moment matrices:

```@docs
ShiftNullspace
Echelon
```

The [`Echelon`](@ref) uses the RowEchelon package to determine the independent
rows (which is not numerically stable) while the [`ShiftNullspace`](@ref) uses
[`RankCheck`](@ref)s with the singular values so it should have better numerical
behavior. They can either simply distinguish the dependency of rows with
[`AnyDependence`](@ref) or use a sieve with [`StaircaseDependence`](@ref) to
save some the computation of the singular values for some submatrices:

```@docs
AnyDependence
StaircaseDependence
```

The relationship between the dependent and the independent rows are
then stored in a [`BorderBasis`](@ref):

```@docs
BorderBasis
```

Once the center of the atoms are determined, a linear system is solved to determine
the weights corresponding to each dirac.
By default, [`MomentMatrixWeightSolver`](@ref) is used by [`atomic_measure`](@ref) so that if there are small differences between moment values corresponding to the same monomial in the matrix
(which can happen if these moments were computed numerically by a semidefinite proramming solvers, e.g., with [SumOfSquares](https://github.com/jump-dev/SumOfSquares.jl)),
the linear system handles that automatically.
```@docs
MomentMatrixWeightSolver
MomentVectorWeightSolver
```
