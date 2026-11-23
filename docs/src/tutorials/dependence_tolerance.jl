# # Tolerance of the numerical dependence of monomials

# Both the extraction of atoms from a moment matrix and the reduction of
# polynomials modulo a truncated ideal end up with the same numerical
# question: given a matrix ``Z`` with orthonormal columns whose rows are
# indexed by monomials, which monomials are linearly dependent on the
# previous ones?
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

# We compare the following candidates for the tolerance where `σ` are the
# singular values, `r` is the rank and `n` is the number of monomials. The
# first one was used by [`Echelon`](@ref) and the second one is the one
# suggested by perturbation theory as we will see below. The last one is the
# one currently used by the package with the [`accuracy`](@ref) of the
# [`MacaulayNullspace`](@ref).

const CANDIDATES = [
    "σ[r+1]/σ[1]/√n" => (σ, r, n, _) -> max(σ[r+1] / σ[1], eps()) / sqrt(n),
    "σ[r+1]/σ[r]·√n" => (σ, r, n, _) -> max(σ[r+1] / σ[r], eps()) * sqrt(n),
    "accuracy·√n" => (σ, r, n, accuracy) -> accuracy * sqrt(n),
]

# For each case, we report the noise, the signal, and whether each candidate
# gives the true pivots.

function row(name, null::MM.MacaulayNullspace, σ, r, piv)
    A = Matrix(transpose(null.matrix))
    n = size(null.matrix, 1)
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
    for (_, tol) in CANDIDATES
        s *= pivots(A, tol(σ, r, n, null.accuracy)) == piv ? " ✓ |" : " ✗ |"
    end
    ## Check that the package uses the last candidate
    d = MM.BasisDependence{MM.LinearDependence}(null)
    @assert Set(findall(!MM.is_dependent, d.dependence)) ==
            pivots(A, last(CANDIDATES)[2](σ, r, n, null.accuracy))
    return s
end

function table(rows)
    header = "| Case | n | σ[r+1]/σ[1] | σ[r+1]/σ[r] | noise | signal |"
    for (name, _) in CANDIDATES
        header *= " $name |"
    end
    sep = "|" * repeat(" --- |", 6 + length(CANDIDATES))
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
    null = MM.MacaulayNullspace(M)
    A = float(Matrix(SparseArrays.sparse(M)))
    σ = svdvals(A)
    σ = [σ; zeros(size(A, 2) - length(σ))]
    r = size(A, 2) - size(null.matrix, 2)
    F = svd(big.(A); full = true)
    piv = true_pivots(Matrix(F.Vt[(r+1):end, :]))
    return row(join(equalities, ", ") * ", D = $D", null, σ, r, piv)
end

table([
    macaulay([x^2 + y^2 - 1.0, x * y - z], 5),
    macaulay([x^2 + 1 / 3 * y^2 - 0.7, x * y - 1 / 7 * z^3], 5),
    macaulay([x^2 + 1e-3 * y^2 - 0.7, x * y - 1e3 * z^3], 6),
    macaulay([x^3 - 0.3x * y + 1 / 7, y^2 - 3x * z, z^2 - x - 0.25], 5),
])

# The noise is at the level of the round-off errors but it is not at the
# level of `σ[r+1] / σ[1]`: it follows `σ[r+1] / σ[r]` instead. This is
# visible for the third case whose equations are badly scaled: the smallest
# nonzero singular value `σ[r]` is small compared to `σ[1]` so the noise is
# about `1e-13` while `σ[r+1] / σ[1]` is about `1e-16`.
#
# This is explained by the perturbation theory of the singular value
# decomposition [Wed72]: if the matrix is perturbed by a matrix of norm
# ``\epsilon``, the sine of the largest principal angle between the space
# spanned by the singular vectors of the first ``r`` singular values and the
# one of the perturbed matrix is at most ``\epsilon / \delta`` where
# ``\delta`` is the gap between ``\sigma_r`` and the next singular value of
# the perturbed matrix. With ``\epsilon = \sigma_{r+1}``, this is about
# ``\sigma_{r+1} / \sigma_r``.
#
# [Wed72] Wedin, Per-Åke.
# *Perturbation bounds in connection with singular value decomposition.*
# BIT Numerical Mathematics 12 (1972): 99-111.
#
# The Gaussian elimination then amplifies these errors. The backward error
# analysis of Gaussian elimination with partial pivoting gives a bound
# proportional to ``n`` times the growth factor [H02, Chapter 9] but, as
# rounding errors are roughly independent, ``\sqrt{n}`` is a more realistic
# estimate [H02, Section 2.8] and the growth factor of partial pivoting is
# small in practice [TS90]. These analyses consider the rounding errors made
# during the elimination while here, the dominant error is the error of the
# space itself, so the factor ``\sqrt{n}`` was in the end chosen empirically,
# see [MultivariateMoments#111](https://github.com/JuliaAlgebra/MultivariateMoments.jl/pull/111).
# It is enough in all the cases of this tutorial.
#
# [H02] Higham, Nicholas J.
# *Accuracy and Stability of Numerical Algorithms.*
# 2nd edition, SIAM, 2002.
#
# [TS90] Trefethen, Lloyd N. and Schreiber, Robert S.
# *Average-case stability of Gaussian elimination.*
# SIAM Journal on Matrix Analysis and Applications 11.3 (1990): 335-360.
#
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
# the atoms are exactly on the variety.

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
    null = MM.MacaulayNullspace(
        chol.L,
        MB.SubBasis{MB.Monomial}(monos),
        MM.accuracy(chol),
    )
    return row(
        @sprintf("%s, δ = %.0e", name, δ),
        null,
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

# Again, the noise follows `σ[r+1] / σ[r]` and not `σ[r+1] / σ[1]`. The
# tolerance `σ[r+1] / σ[1] / √n` previously used by [`Echelon`](@ref) is
# therefore too small: it considers some noise as pivots. Note that this only
# matters if a monomial is dependent before the last pivot. This is the case
# for atoms on the circle as ``x^2`` is dependent on ``1`` and ``y^2``.
# The last row is discussed in the limitations below.

# ## Moment matrices without gap
#
# Should we then simply use `σ[r+1] / σ[r]` ? Consider the following moment
# matrix obtained by an SDP solver for [LPJ20, Example 3.9].
#
# [LPJ20] Legat, Benoît, Pablo Parrilo, and Raphaël Jungers.
# *Certifying unstability of switched systems using sum of squares
# programming.*
# SIAM Journal on Control and Optimization 58.4 (2020): 2616-2638.

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
σ = svdvals(Matrix(MM.value_matrix(ν)))
σ / σ[1]

# The singular values decrease without any gap. Whatever the rank we choose,
# the discarded singular values are not noise, they are just smaller. The
# space spanned by the first `r` singular vectors is then very sensitive to
# perturbations and `σ[r+1] / σ[r]` is large. Still, the atoms can be
# extracted as long as the Gaussian elimination does not consider the small
# but nonzero pivots as zero. We compare the atoms obtained with the
# tolerance `σ[r+1] / σ[r] · √n` and with the one of the package.

function extract(rank_tol, accuracy)
    null = MM.image_space(ν, MM.LeadingRelativeRankTol(rank_tol))
    r = size(null.matrix, 2)
    if !isnothing(accuracy)
        null = MM.MacaulayNullspace(null.matrix, null.basis, accuracy(r))
    end
    atoms = MM.solve(null, MM.Echelon())
    if isnothing(atoms) || isempty(atoms)
        return r, "none"
    end
    return r,
    join([@sprintf("%.4f", atom[2] / atom[1]) for atom in atoms], ", ")
end

rows = map([1e-4, 1e-5, 1e-6, 1e-7, 1e-8]) do rank_tol
    r, gap = extract(rank_tol, r -> σ[r+1] / σ[r])
    _, package = extract(rank_tol, nothing)
    return @sprintf(
        "| %.0e | %d | %.1e | %.1e | %s | %s |",
        rank_tol,
        r,
        σ[r+1] / σ[1],
        σ[r+1] / σ[r],
        gap,
        package,
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

# The dominant atom has `b / a ≈ 0.2278`. With `σ[r+1] / σ[r]`, the tolerance
# is so large that it considers true pivots as zero and the extraction fails.

# ## The first gap
#
# The two previous sections seem contradictory. The difference is whether the
# singular values discarded by the rank check are noise:
#
# * For a Macaulay matrix of exact equations or a moment matrix with a clear
#   rank, there is a gap right after `σ[r]` and `σ[r+1]` is the level of the
#   noise. The error of the space is then `σ[r+1] / σ[r]`.
# * For the moment matrix above, there is no gap so the level of the noise
#   cannot be determined from the singular values. The space spanned by the
#   first `r` singular vectors should then be considered as the object of
#   interest rather than an approximation of the space of a matrix of rank
#   `r`. The Gaussian elimination should only account for the errors in
#   computing this space, not for the distance to the space of a matrix of
#   rank `r`, hence the use of `σ[r+1] / σ[1]`.
#
# These are two extremes of the same rule. Let `j ≥ r` be the first index
# such that there is a gap after `σ[j]`, i.e., `σ[j+1] < 10⁻² σ[j]`. The
# singular values `σ[r+1], ..., σ[j]` were discarded by the rank check but they
# are not noise, the level of the noise is `σ[j+1]`. The error of the space
# is then `σ[j+1] / σ[r]`. If there is no gap, we fall back to
# `σ[r+1] / σ[1]`. This is what [`accuracy`](@ref) computes.
#
# The fallback is also useful when the rank is chosen larger than the rank
# supported by the data. Consider 6 atoms on the twisted cubic, two of which
# are very close to each other so that `σ[6]` is about `1e-10 σ[1]`.

twisted_cubic_atoms =
    [[t, t^2, t^3] for t in big.([-0.9, -0.4, 0.1, 0.5, 0.5001, 0.95])]
table([
    moments("twisted cubic, d = 3", [x, y, z], 3, twisted_cubic_atoms, δ) for
    δ in [1e-12, 1e-8, 1e-6]
])

# The first 6 monomials are independent so there is no noise: only the signal
# matters. With `δ = 1e-8` and `δ = 1e-6`, `σ[6]` is at the level of the noise
# so `σ[7] / σ[6]` is close to `1` and the tolerance `σ[r+1] / σ[r] · √n`
# considers true pivots as zero. As there is no gap after `σ[6]`, the
# accuracy falls back to `σ[7] / σ[1]` which gives the true pivots. With
# `δ = 1e-12`, `σ[7] / σ[6]` is just above the threshold of `10⁻²` so the
# accuracy also falls back to `σ[7] / σ[1]`.
#
# ## Limitations
#
# The threshold `10⁻²` for detecting a gap is arbitrary. When the noise of an
# SDP solution is large, for instance `δ = 1e-4` for the sphere above,
# `σ[r+1] / σ[r]` is about `0.2` which is not considered as a gap so the
# accuracy falls back to `σ[r+1] / σ[1]` which is too small. Increasing the
# threshold to detect it would also consider `σ[3] / σ[2] ≈ 0.096` of the
# moment matrix of the previous section as a gap. Even then, the tolerance
# `σ[r+1] / σ[r] · √n` would be larger than the signal. The true pivots are
# obtained with a tolerance between the noise `1.3e-2` and the signal
# `1.1e-1` but the singular values alone do not tell us where these are. In
# general, when the noise is comparable to the signal, the dependence of the
# monomials cannot be recovered with any tolerance.
