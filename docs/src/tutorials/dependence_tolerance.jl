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

# For each case, the first table gives the number of monomials `n`, the
# measured noise and signal, and two ratios. With `r` the selected rank,
# `ρ₁ = σ[r+1] / σ[1]` scales the first discarded singular value by the
# largest one, and `ρᵣ = σ[r+1] / σ[r]` scales it by the weakest retained one.
#
# The recovery summary counts exact matches with the reference pivot set.
# Methods share a row when they succeed and fail on the same cases in that
# experiment. Two methods use fixed tolerances, `ρ₁/√n` (previously used by
# [`Echelon`](@ref)) and `ρᵣ·√n`. Both floor the ratio at machine precision
# before scaling and apply a single elimination. **Gap at r** checks for
# a gap at the selected rank, and **First gap** uses the first gap at or
# after that rank. Both call [`BasisDependence`](@ref) and retry at a
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

function row(labels, Z, basis, σ, r, piv)
    A = Matrix(transpose(Z))
    n = size(Z, 1)
    noise, signal = noise_and_signal(A, piv)
    label = join(labels, " | ")
    measurements = @sprintf(
        "| %s | %d | %.1e | %.1e | %.1e | %.1e |",
        label,
        n,
        σ[r+1] / σ[1],
        σ[r+1] / σ[r],
        noise,
        signal,
    )
    leading = max(σ[r+1] / σ[1], eps()) / sqrt(n)
    adjacent = max(σ[r+1] / σ[r], eps()) * sqrt(n)
    at_rank, first_gap = gap_pivots(Z, basis, σ, r)
    found = (pivots(A, leading), pivots(A, adjacent), at_rank, first_gap)
    return (; measurements, case = join(labels, " "), reference = piv, found)
end

function recovery_summary(rows)
    methods = ["ρ₁/√n", "ρᵣ·√n", "Gap at r", "First gap"]
    cases = [row.case for row in rows]
    outcomes = [
        [row.found[i] == row.reference for row in rows] for
        i in eachindex(methods)
    ]
    return map(unique(outcomes)) do pattern
        names = join(methods[findall(==(pattern), outcomes)], ", ")
        missed = if all(pattern)
            "None"
        elseif any(pattern)
            join(cases[.!pattern], ", ")
        else
            "All"
        end
        return "| $names | $(count(pattern))/$(length(rows)) | $missed |"
    end
end

function table(rows, labels = ("Case",))
    header = join(labels, " | ")
    sep = repeat(" --- |", length(labels))
    return Markdown.parse("""
    | $header | n | ρ₁ | ρᵣ | Noise | Signal |
    |$sep --- | --- | --- | --- | --- |
    $(join((row.measurements for row in rows), "\n"))

    **Recovery of the reference pivots**

    | Method(s) | Recovered | Missed cases |
    | --- | --- | --- |
    $(join(recovery_summary(rows), "\n"))
    """)
end

# ## Macaulay null spaces
#
# The Macaulay matrix of degree `D` of exact equations has entries known up
# to machine precision so its singular values have a clear gap between the
# nonzero singular values and the ones that are zero up to round-off errors.

@polyvar x y z

function macaulay(label, equalities, D)
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
    return row((label,), Z, M.column_basis, F.S, r, piv)
end

# The four systems below are truncated at degree `D`.
#
# * **M1** (D = 5) uses `x² + y² = 1` and `xy = z`.
# * **M2** (D = 5) uses `x² + y²/3 = 0.7` and `xy = z³/7`.
# * **M3** (D = 6) uses `x² + 0.001y² = 0.7` and `xy = 1000z³`.
# * **M4** (D = 5) uses `x³ − 0.3xy + 1/7 = 0`, `y² = 3xz`, and
#   `z² = x + 0.25`.

table([
    macaulay("M1", [x^2 + y^2 - 1.0, x * y - z], 5),
    macaulay("M2", [x^2 + 1 / 3 * y^2 - 0.7, x * y - 1 / 7 * z^3], 5),
    macaulay("M3", [x^2 + 1e-3 * y^2 - 0.7, x * y - 1e3 * z^3], 6),
    macaulay("M4", [x^3 - 0.3x * y + 1 / 7, y^2 - 3x * z, z^2 - x - 0.25], 5),
])

# In M3 the equations have very different scales. The noise
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
# an SDP solver, we perturb the moment matrix by a random symmetric matrix.
# For each set of atoms we hold the perturbation direction fixed and vary
# only its relative norm `δ`. The rank is fixed to the number of atoms with
# [`FixedRank`](@ref) to separate the choice of the rank from the choice of
# the pivots. The atoms and the true pivots are computed in `BigFloat` so that
# roundoff in the reference pivots is small compared with the perturbation.

function moments(vars, d, atoms, levels; case = nothing)
    monos = monomials(vars, 0:d)
    k = length(atoms)
    V = [m(vars => a) for m in monos, a in atoms]
    W = big.(range(1, 2; length = k))
    exact = Float64.(V * Diagonal(W) * V')
    piv = true_pivots(Matrix(transpose(V)))
    E = randn(size(exact))
    E = (E + E') / 2
    scale = opnorm(exact) / opnorm(E)
    basis = MB.SubBasis{MB.Monomial}(monos)
    return map(levels) do δ
        M = exact + (δ * scale) * E
        chol = MM.low_rank_ldlt(M, MM.SVDLDLT(), MM.FixedRank(k))
        noise_level = @sprintf("%.0e", δ)
        labels = isnothing(case) ? (noise_level,) : (case, noise_level)
        return row(labels, chol.L, basis, chol.singular_values, k, piv)
    end
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

# All three cases use monomials through degree four. **C₁** and **C₃** each
# have seven atoms on a circle of radius one and three, respectively.
# **S** has twelve atoms on the unit sphere. Each row gives the relative
# perturbation `δ` for that case.

Random.seed!(0)
circle_atoms = circle(1, 7)
large_circle_atoms = circle(3, 7)
sphere_atoms = [sphere() for _ in 1:12]
table(
    vcat(
        moments(
            [x, y],
            4,
            circle_atoms,
            [1e-12, 1e-8, 1e-6, 1e-4];
            case = "C₁",
        ),
        moments(
            [x, y],
            4,
            large_circle_atoms,
            [1e-12, 1e-8, 1e-6];
            case = "C₃",
        ),
        moments(
            [x, y, z],
            4,
            sphere_atoms,
            [1e-12, 1e-8, 1e-6, 5e-6, 6e-6, 1e-4];
            case = "S",
        ),
    ),
    ("Case", "δ"),
)

# On both circles there is already a gap at the retained rank. Both gap
# policies therefore start at `ρᵣ·√n` (the roundoff floor is inactive).
# The smaller `ρ₁/√n` admits noise in dependent columns such as `x²`,
# which on the circle is a linear combination of `1` and `y²`.
#
# The two nearby sphere perturbations show what happens at the gap cutoff.
# At `δ = 5e-6`, `ρᵣ ≈ 0.0099` qualifies as a gap. All three methods using
# `ρᵣ·√n` recover the reference pivots. At `δ = 6e-6`, `ρᵣ ≈ 0.0118`
# exceeds the cutoff and there is no qualifying gap later. The fixed
# tolerance is about `0.070`, between the measured noise (`0.0020`) and
# signal (`0.114`). The two gap policies fall back to `ρ₁·√n ≈ 0.000029`,
# accept column 10 in place of column 13, and still find all twelve pivots.
#
# At `δ = 1e-4` all four methods miss the reference. The fixed adjacent
# tolerance exceeds one, so it discards every entry of the orthonormal
# basis. Both gap policies start below the measured noise and again accept
# a spurious pivot. [Choosing which gap to use](@ref) compares matrices
# where the choice of search range matters.

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
# and when we use the package default. The last two columns give the atom
# ratios `b/a` for these two choices. Both use the same computed vectors
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
            "| Rank tol. | r | ρ₁ | ρᵣ | b/a (ρᵣ) | b/a (default) |"
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
# **A. Spurious coordinate.** In the diagonal matrix below, the clean
# rank-two subspace is spanned by coordinate vectors 1 and 3. A small
# perturbation couples directions 2 and 3 while leaving the fourth singular
# value zero.

A = diagm([1.0, 5e-4, 1e-3, 0.0])
A[2, 3] = A[3, 2] = 1e-12
F_coupled = svd(A)
Z_coupled = F_coupled.U[:, 1:2]
basis_coupled = MB.SubBasis{MB.Monomial}(monomials([x], 0:3))
coupled = row(("A",), Z_coupled, basis_coupled, F_coupled.S, 2, Set([1, 3]))
@assert coupled.found[3:4] == (Set([1, 3]), Set([1, 2]))

# The selected-rank policy skips the spurious `2e-9` component in column 2.
# Searching for a later gap finds one after the third singular value and
# sets a tolerance of about `1.8e-12`, which accepts that component.
#
# **B. Genuine coordinate.** The following moment matrix comes from
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
two_atoms = row(("B",), Z_atoms, basis_atoms, F_atoms.S, 1, Set([1]))
@assert two_atoms.found[3:4] == (Set([4]), Set([1]))

# Each entry gives the column indices selected in the displayed basis.

function pivot_row(result)
    cells = map((result.reference, result.found...)) do indices
        return isempty(indices) ? "∅" : join(sort(collect(indices)), ", ")
    end
    return "| $(result.case) | $(join(cells, " | ")) |"
end

Markdown.parse("""
| Case | Reference | ρ₁/√n | ρᵣ·√n | Gap at r | First gap |
| --- | --- | --- | --- | --- | --- |
$(pivot_row(coupled))
$(pivot_row(two_atoms))
""")

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
# We use monomials through degree three and vary the relative perturbation `δ`.

twisted_cubic_atoms =
    [[t, t^2, t^3] for t in big.([-0.9, -0.4, 0.1, 0.5, 0.5001, 0.95])]
table(moments([x, y, z], 3, twisted_cubic_atoms, [1e-12, 1e-8, 1e-6]), ("δ",))

# The first six monomials are already independent. A generic six-dimensional
# subspace has the same six pivots, so matching the reference says little about
# recovery of the clean subspace. At higher perturbations the sixth
# singular direction can be poorly determined.
#
# The gap threshold `10⁻²` is heuristic. For the sphere at `δ = 1e-4`, the
# fallback can accept spurious pivots. The noise and signal columns show
# the interval of tolerances that recovers the reference pivots. Once the
# measured noise exceeds the signal, no single tolerance along this
# elimination path can recover them all.
