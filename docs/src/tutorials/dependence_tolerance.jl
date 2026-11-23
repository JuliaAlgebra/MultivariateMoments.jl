# # Tolerance of the numerical dependence of monomials

# In the SVD examples below, atom extraction from a moment matrix and
# reduction modulo a truncated ideal give a matrix ``Z`` with orthonormal
# columns. Its rows are indexed by monomials. Which rows are linearly
# dependent on the previous ones?
#
# * For atom extraction, ``Z`` spans the image space of the moment matrix,
#   see `image_space`. The moment matrix is typically the solution of
#   an SDP solver hence it is only known with the accuracy of the solver.
# * For a truncated ideal, ``Z`` spans the null space of the Macaulay matrix
#   of exact equations, see [`MacaulayNullspace`](@ref). It is used,
#   e.g., by SumOfSquares to reduce the Gram basis of a Sum-of-Squares
#   certificate over an algebraic set.
#
# In both cases, ``Z`` is computed from a singular value decomposition and
# the dependent monomials are the columns that are not pivots of the reduced
# row echelon form of ``Z^\top``. The Gaussian elimination needs a tolerance
# below which an entry is considered to be zero. This tutorial studies how to
# choose this tolerance so that it works for both cases. It explains the
# choice made in [`accuracy`](@ref).

using LinearAlgebra, Markdown, Printf, Random
import GenericLinearAlgebra # `svd` of `Matrix{BigFloat}`
using DynamicPolynomials
import MultivariateBases as MB
import MultivariateMoments as MM
import SparseArrays

# ## Noise and signal
#
# We compute the true pivots by doing the same computation in `BigFloat`.
# To make sure the reference is not affected by the round-off errors, the
# elimination in `BigFloat` uses a tolerance much smaller than the precision
# of `Float64` but much larger than the precision of `BigFloat`.

setprecision(BigFloat, 128)

function pivots(A, tol)
    U = copy(A)
    MM.RowEchelon.rref!(U, tol)
    return Set(MM._pivots(U))
end

true_pivots(A::Matrix{BigFloat}) = pivots(A, big(1e-25))

# Given the true pivots, we run the Gaussian elimination by following them
# instead of comparing the entries with a tolerance. Along the way, we record:
#
# * the *noise*: the largest entry of a column that is truly dependent; it
#   should be zero but it is nonzero because of numerical errors,
# * the *signal*: the smallest pivot of a column that is truly independent.
#
# A tolerance `tol` gives the true pivots if `noise < tol < signal`.

function noise_and_signal(A, piv)
    A = copy(A)
    nr, nc = size(A)
    i = 1
    noise = 0.0
    signal = Inf
    for j in 1:nc
        if i > nr
            break
        end
        m, mi = findmax(abs, view(A, i:nr, j))
        mi += i - 1
        if j in piv
            signal = min(signal, m)
            A[[i, mi], :] = A[[mi, i], :]
            A[i, :] ./= A[i, j]
            for k in 1:nr
                if k != i
                    A[k, :] .-= A[k, j] .* A[i, :]
                end
            end
            i += 1
        else
            noise = max(noise, m)
            A[i:nr, j] .= 0
        end
    end
    return noise, signal
end

# We compare two fixed tolerances with the pivots selected by the package.
# Here `σ` contains the singular values, `r` is the selected rank, and `n` is
# the number of monomials. The first tolerance was previously used by
# [`Echelon`](@ref). The second uses the smallest retained singular value to
# scale the first discarded one. Both apply a single Gaussian elimination.
# The last two columns use the gap at the selected rank or search for a
# later gap, respectively. Both call [`BasisDependence`](@ref) and retry at a
# lower tolerance if the initial elimination loses rank. For example,
# a rank-one space spanned by a vector with 128 equal entries has unit
# norm after normalization, but each entry is only about 0.088. If the
# proposed tolerance is 0.09, the first pass discards every entry.

function gap_pivots(Z, basis, σ, r; gap_tol = 1e-2)
    return map((false, true)) do search_later
        policy = MM.GapBasedAccuracy(MM.FixedRank(r); gap_tol, search_later)
        null = MM.MacaulayNullspace(Z, basis, MM.accuracy(σ, r, policy))
        d = MM.BasisDependence{MM.LinearDependence}(null)
        return Set(findall(!MM.is_dependent, d.dependence))
    end
end

function row(name, Z, basis, σ, r, piv)
    A = Matrix(transpose(Z))
    n = size(Z, 1)
    noise, signal = noise_and_signal(A, piv)
    s = @sprintf(
        "| %s | %d | %.1e | %.1e | %.1e | %.1e |",
        name,
        n,
        σ[r+1] / σ[1],
        σ[r+1] / σ[r],
        noise,
        signal,
    )
    for tol in (
        max(σ[r+1] / σ[1], eps()) / sqrt(n),
        max(σ[r+1] / σ[r], eps()) * sqrt(n),
    )
        s *= pivots(A, tol) == piv ? " ✓ |" : " ✗ |"
    end
    for found in gap_pivots(Z, basis, σ, r)
        s *= found == piv ? " ✓ |" : " ✗ |"
    end
    return s
end

function table(rows)
    header = "| Case | n | σ[r+1]/σ[1] | σ[r+1]/σ[r] | noise | signal | σ[r+1]/σ[1]/√n | σ[r+1]/σ[r]·√n | Gap at r | First gap at or after r |"
    sep = "|" * repeat(" --- |", 10)
    return Markdown.parse(join([header; sep; rows], "\n"))
end

# ## Macaulay null spaces
#
# The Macaulay matrix of degree `D` of exact equations has entries known up
# to machine precision so its singular values have a clear gap between the
# nonzero singular values and the ones that are zero up to round-off errors.

@polyvar x y z

function macaulay(equalities, D)
    M = MB.macaulay(equalities, D, [x, y, z])
    A = float(Matrix(SparseArrays.sparse(M)))
    F = svd(A; full = true)
    check = MM._default_nullspace_rank_check(A)
    r = MM.rank_from_singular_values(F.S, check)
    Z = F.Vt[(r+1):end, :]'
    big_A = big.(A)
    ref = svd(big_A; full = true)
    ref_rank = MM.rank_from_singular_values(
        ref.S,
        MM._default_nullspace_rank_check(big_A),
    )
    @assert r == ref_rank
    piv = true_pivots(Matrix(ref.Vt[(ref_rank+1):end, :]))
    return row(
        join(equalities, ", ") * ", D = $D",
        Z,
        M.column_basis,
        F.S,
        r,
        piv,
    )
end

table([
    macaulay([x^2 + y^2 - 1.0, x * y - z], 5),
    macaulay([x^2 + 1 / 3 * y^2 - 0.7, x * y - 1 / 7 * z^3], 5),
    macaulay([x^2 + 1e-3 * y^2 - 0.7, x * y - 1e3 * z^3], 6),
    macaulay([x^3 - 0.3x * y + 1 / 7, y^2 - 3x * z, z^2 - x - 0.25], 5),
])

# In the third example the equations have very different scales. The noise
# in the elimination is about `1e-13` even though `σ[r+1] / σ[1]` is close to
# machine precision. Scaling by `σ[r]` gives a larger tolerance.
#
# Wedin bounds the angle between singular subspaces using perturbation or
# residual norms and the spectral separation of retained and discarded
# directions. With a pronounced gap at rank `r`, that separation may be
# close to `σ[r]`. Using `σ[r+1]` as a proxy for the perturbation size then
# suggests the scale `σ[r+1] / σ[r]`. An actual angle bound needs independent
# information about the perturbation. [Wedin1972](@cite)
#
# Gaussian elimination introduces rounding error and can amplify errors
# already present in the singular vectors. For partial pivoting, a worst-case
# backward error bound includes a dimension factor and element growth
# (Chapter 9). The square-root rule of thumb for accumulated rounding error
# assumes cancellation of errors of varying signs (Section 2.8).
# [Higham2002](@cite)
#
# Average-case experiments show much smaller element growth than the
# worst-case bound for several random matrix ensembles.
# [Trefethen1990](@cite)
#
# Applying these results to dependence decisions requires an additional
# empirical choice. Here the input may already contain error from an SDP
# solver. The multiplier `√n` for the initial pivot tolerance was chosen
# empirically in [pull request #111](https://github.com/JuliaAlgebra/MultivariateMoments.jl/pull/111).
# The examples below compare the two choices of spectral gap.

# For exact equations, the computed `σ[r+1]` is often below the machine
# precision. It is not a reliable estimate of the noise since the singular
# values are only computed with an error of about `length(σ) * eps * σ[1]`,
# which is the default tolerance of `LinearAlgebra.rank`. For instance, for
# the equations `x^2 + y^2 / 3 - 0.7` and `x * y - y^3 / 7` at degree 6, the
# computed `σ[r+1] / σ[1]` changes from `1.0e-16` to `7.4e-17` between Julia
# v1.13 and v1.15, and the noise of the Gaussian elimination, about
# `1.9e-15`, exceeds `σ[r+1] / σ[r] · √n` with v1.15. The level of noise used
# by [`accuracy`](@ref) is therefore at least `length(σ) * eps * σ[1]`.

# ## Moment matrices with noise
#
# We now consider the moment matrix of atoms lying on a variety so that some
# monomials are dependent before the last pivot. To mimic the inaccuracy of
# an SDP solver, the moment matrix is perturbed by a random symmetric matrix of
# relative norm `δ`. The rank is fixed to the number of atoms with
# [`FixedRank`](@ref) to separate the choice of the rank from the choice of
# the pivots. The atoms and the true pivots are computed in `BigFloat` so that
# roundoff in the reference pivots is small compared with the perturbation.

function moments(name, vars, d, atoms, δ)
    monos = monomials(vars, 0:d)
    k = length(atoms)
    V = [m(vars => a) for m in monos, a in atoms]
    W = big.(range(1, 2; length = k))
    exact = Float64.(V * Diagonal(W) * V')
    piv = true_pivots(Matrix(transpose(V)))
    E = randn(size(exact))
    E = (E + E') / 2
    M = exact + δ * opnorm(exact) / opnorm(E) * E
    chol = MM.low_rank_ldlt(M, MM.SVDLDLT(), MM.FixedRank(k))
    return row(
        @sprintf("%s, δ = %.0e", name, δ),
        chol.L,
        MB.SubBasis{MB.Monomial}(monos),
        chol.singular_values,
        k,
        piv,
    )
end

## Atoms on the circle of radius `radius` with angles `2π i / k + 0.1 i^2` for
## `i = 1, ..., k`
function circle(radius, k)
    return [
        radius .* [cos(t), sin(t)] for
        t in (2big(π) * i / k + big(0.1) * i^2 for i in 1:k)
    ]
end
function sphere()
    u, v = big(π) * rand() - big(π) / 2, 2big(π) * rand()
    return [cos(u) * cos(v), cos(u) * sin(v), sin(u)]
end

Random.seed!(0)
circle_atoms = circle(1, 7)
large_circle_atoms = circle(3, 7)
sphere_atoms = [sphere() for _ in 1:12]
table(
    [
        [
            moments("circle, d = 4", [x, y], 4, circle_atoms, δ) for
            δ in [1e-12, 1e-8, 1e-6, 1e-4]
        ]
        [
            moments(
                "circle of radius 3, d = 4",
                [x, y],
                4,
                large_circle_atoms,
                δ,
            ) for δ in [1e-12, 1e-8, 1e-6]
        ]
        [
            moments("sphere, d = 4", [x, y, z], 4, sphere_atoms, δ) for
            δ in [1e-12, 1e-8, 1e-6, 1e-4]
        ]
    ],
)

# On the circle, `x²` is dependent on `1` and `y²`. A tolerance below the
# measured noise can select this column as a pivot. Several rows illustrate
# this failure for the first tolerance. The last row is discussed below.

# ## Moment matrices without gap
#
# Should we then simply use `σ[r+1] / σ[r]` ? Consider the following moment
# matrix obtained by an SDP solver for [Legat2020; Example 3.9](@cite).

@polyvar a b
ν = MM.moment_matrix(
    [
        14.8107 3.37426 0.769196 0.175447 0.0400656 0.00917416 0.00211654
        3.37426 0.769196 0.175447 0.0400656 0.00917416 0.00211654 0.000498924
        0.769196 0.175447 0.0400656 0.00917416 0.00211654 0.000498924 0.000125543
        0.175447 0.0400656 0.00917416 0.00211654 0.000498924 0.000125543 3.76638e-5
        0.0400656 0.00917416 0.00211654 0.000498924 0.000125543 3.76638e-5 1.62883e-5
        0.00917416 0.00211654 0.000498924 0.000125543 3.76638e-5 1.62883e-5 1.07817e-5
        0.00211654 0.000498924 0.000125543 3.76638e-5 1.62883e-5 1.07817e-5 1.10658e-5
    ][
        7:-1:1,
        7:-1:1,
    ],
    monomials([a, b], 6),
)
F = svd(Matrix(MM.value_matrix(ν)))
σ = F.S
σ / σ[1]

# After the first singular value the spectrum decays gradually. The table
# compares the extracted atoms when we set the accuracy to `σ[r+1] / σ[r]`
# and when we use the package default. Both use the same computed vectors
# and the [`Echelon`](@ref) solver.

function extract(null)
    atoms = MM.solve(null, MM.Echelon())
    if isnothing(atoms) || isempty(atoms)
        return "none"
    end
    return join([@sprintf("%.4f", atom[2] / atom[1]) for atom in atoms], ", ")
end

rows = map([1e-4, 1e-5, 1e-6, 1e-7, 1e-8]) do rank_tol
    check = MM.LeadingRelativeRankTol(rank_tol)
    r = MM.rank_from_singular_values(σ, check)
    Z = F.U[:, 1:r]
    adjacent = MM.MacaulayNullspace(Z, ν.basis, σ[r+1] / σ[r])
    package = MM.MacaulayNullspace(Z, ν.basis, MM.accuracy(σ, r, check))
    return @sprintf(
        "| %.0e | %d | %.1e | %.1e | %s | %s |",
        rank_tol,
        r,
        σ[r+1] / σ[1],
        σ[r+1] / σ[r],
        extract(adjacent),
        extract(package),
    )
end

Markdown.parse(
    join(
        [
            "| Rank tolerance | r | σ[r+1]/σ[1] | σ[r+1]/σ[r] | b/a with σ[r+1]/σ[r] | b/a with `accuracy` |"
            "| --- | --- | --- | --- | --- | --- |"
            rows
        ],
        "\n",
    ),
)

# The dominant atom has `b / a ≈ 0.2278`. Lowering the rank tolerance retains
# additional directions and makes the choice of pivots relevant to which
# atoms can be extracted.

# ## Choosing which gap to use
#
# [`GapBasedAccuracy`](@ref) wraps a rank check. With
# `gap_tol = 1e-2`, a gap after `j` means `σ[j+1] < 10⁻² σ[j]`.
# The default, `search_later = false`, checks only `j = r`. Setting
# `search_later = true` uses the first qualifying gap at or after `r`. If
# either policy finds a gap, it divides `σ[j+1]` by `σ[r]`. Otherwise it
# divides `σ[r+1]` by `σ[1]`. The numerator includes the roundoff floor
# described above. A later gap can suggest a noise floor for a larger signal
# space. Wedin's separation for the chosen subspace is still the separation
# at rank `r`. Changing the search range can even increase the tolerance.
# For example, with `σ = [1, 1e-3, 5e-4, 4e-6]` and `r = 2`, the default
# gives `5e-4`. Searching one step further gives `4e-6 / 1e-3 = 4e-3`.
#
# **A small spurious coordinate.** In the diagonal matrix below, the clean
# rank-two subspace is spanned by coordinate vectors 1 and 3. A small
# perturbation couples directions 2 and 3 while leaving the fourth singular
# value zero.

A = diagm([1.0, 5e-4, 1e-3, 0.0])
A[2, 3] = A[3, 2] = 1e-12
F_coupled = svd(A)
Z_coupled = F_coupled.U[:, 1:2]
basis_coupled = MB.SubBasis{MB.Monomial}(monomials([x], 0:3))
coupled_pivots = gap_pivots(Z_coupled, basis_coupled, F_coupled.S, 2)
@assert coupled_pivots == (Set([1, 3]), Set([1, 2]))

# The selected-rank policy skips the spurious `2e-9` component in column 2.
# Searching for a later gap finds one after the third singular value and
# sets a tolerance of about `1.8e-12`, which accepts that component.
#
# **A small genuine coordinate.** The following moment matrix comes from
# the package's atom-extraction tests. Its measure has atoms at `(1, 3)`
# and `(4, 3)` with approximate weights `0.864` and `0.136`. We retain one
# direction of this two-atom example. The clean moment matrix has strictly
# positive entries, so its dominant eigenvector can be chosen strictly
# positive. Its first coordinate, the constant monomial, is a genuine pivot.

Q = [
    1.0000000000000002 2.9999999777389235 1.408602019734018 8.999999911981313 4.225805987406138 3.043010098670093
    2.9999999777389235 8.999999911981313 4.225805987406138 26.999999824988187 12.677417999642149 9.129030026074997
    1.408602019734018 4.225805987406138 3.043010098670093 12.677417999642149 9.129030026074997 9.580642414414385
    8.999999911981313 26.999999824988187 12.677417999642149 80.99999955990646 38.03225428611012 27.387090350285558
    4.225805987406138 12.677417999642149 9.129030026074997 38.03225428611012 27.387090350285558 28.74192618075039
    3.043010098670093 9.129030026074997 9.580642414414385 27.387090350285558 28.74192618075039 35.73117167739159
]
F_atoms = svd(Q)
Z_atoms = F_atoms.U[:, 1:1]
basis_atoms = MB.SubBasis{MB.Monomial}(monomials([x, y], 0:2))
atom_pivots = gap_pivots(Z_atoms, basis_atoms, F_atoms.S, 1)
@assert atom_pivots == (Set([4]), Set([1]))

table([
    row(
        "coupled diagonal matrix",
        Z_coupled,
        basis_coupled,
        F_coupled.S,
        2,
        Set([1, 3]),
    ),
    row("two-atom moment matrix", Z_atoms, basis_atoms, F_atoms.S, 1, Set([1])),
])

# The second singular value is about `0.21 σ[1]`, well above the residual
# singular values from the third onward. The selected-rank fallback sets
# a tolerance of about `0.51`, skipping the genuine first coordinate of
# magnitude `0.083`. The gap after the second direction gives a tolerance
# of about `9.5e-9` and recovers the constant monomial. To use this policy
# when computing the image space, pass it as the rank check.

measure_moments = MM.moment_matrix(Q, monomials([x, y], 0:2))
policy = MM.GapBasedAccuracy(MM.FixedRank(1); search_later = true)
image = MM.image_space(measure_moments, policy)
d = MM.BasisDependence{MM.LinearDependence}(image)
@assert findall(!MM.is_dependent, d.dependence) == [1]

# The cutoff controls which gaps qualify. Increasing it to `0.25` makes
# the gap after the first singular value qualify for both policies.

@assert gap_pivots(Z_atoms, basis_atoms, F_atoms.S, 1; gap_tol = 0.25) ==
        (Set([4]), Set([4]))

# A later gap can help when discarded directions still carry signal, as in
# the two-atom example. A gap far below the selected rank can also suggest
# a tolerance small enough to admit roundoff or perturbation in a dependent
# column, as in the coupled matrix. The clean pivots distinguish these two
# outcomes even though both policies recover the requested number of pivots.

# ## Limitations
#
# Consider six atoms on the twisted cubic, two of which are very close.

twisted_cubic_atoms =
    [[t, t^2, t^3] for t in big.([-0.9, -0.4, 0.1, 0.5, 0.5001, 0.95])]
table([
    moments("twisted cubic, d = 3", [x, y, z], 3, twisted_cubic_atoms, δ) for
    δ in [1e-12, 1e-8, 1e-6]
])

# The first six monomials are already independent. A generic six-dimensional
# subspace has the same six pivots, so a check mark here says little about
# recovery of the clean subspace. At higher perturbations the sixth
# singular direction can be poorly determined.
#
# The gap threshold `10⁻²` is heuristic. For the sphere at `δ = 1e-4`, the
# fallback can accept spurious pivots. The noise and signal columns show
# the interval of tolerances that recovers the reference pivots. Once the
# measured noise exceeds the signal, no single tolerance along this
# elimination path can recover them all.
