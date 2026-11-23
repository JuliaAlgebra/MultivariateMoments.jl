using LinearAlgebra, Markdown, Printf, Random
import GenericLinearAlgebra # `svd` of `Matrix{BigFloat}`
using DynamicPolynomials
import MultivariateBases as MB
import MultivariateMoments as MM
import SparseArrays

setprecision(BigFloat, 128)

function pivots(A, tol)
    U = copy(A)
    MM.RowEchelon.rref!(U, tol)
    return Set(MM._pivots(U))
end

true_pivots(A::Matrix{BigFloat}) = pivots(A, big(1e-25))

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

table([
    macaulay("M1", [x^2 + y^2 - 1.0, x * y - z], 5),
    macaulay("M2", [x^2 + 1 / 3 * y^2 - 0.7, x * y - 1 / 7 * z^3], 5),
    macaulay("M3", [x^2 + 1e-3 * y^2 - 0.7, x * y - 1e3 * z^3], 6),
    macaulay("M4", [x^3 - 0.3x * y + 1 / 7, y^2 - 3x * z, z^2 - x - 0.25], 5),
])

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

# Atoms on the circle of radius `radius` with angles `2π i / k + 0.1 i^2` for
# `i = 1, ..., k`
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

A = diagm([1.0, 5e-4, 1e-3, 0.0])
A[2, 3] = A[3, 2] = 1e-12
F_coupled = svd(A)
Z_coupled = F_coupled.U[:, 1:2]
basis_coupled = MB.SubBasis{MB.Monomial}(monomials([x], 0:3))
coupled = row(("A",), Z_coupled, basis_coupled, F_coupled.S, 2, Set([1, 3]))
@assert coupled.found[3:4] == (Set([1, 3]), Set([1, 2]))

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

measure_moments = MM.moment_matrix(Q, monomials([x, y], 0:2))
policy = MM.GapBasedAccuracy(MM.FixedRank(1); search_later = true)
image = MM.image_space(measure_moments, policy)
d = MM.BasisDependence{MM.LinearDependence}(image)
@assert findall(!MM.is_dependent, d.dependence) == [1]

@assert gap_pivots(Z_atoms, basis_atoms, F_atoms.S, 1; gap_tol = 0.25) ==
        (Set([4]), Set([4]))

twisted_cubic_atoms =
    [[t, t^2, t^3] for t in big.([-0.9, -0.4, 0.1, 0.5, 0.5001, 0.95])]
table(moments([x, y, z], 3, twisted_cubic_atoms, [1e-12, 1e-8, 1e-6]), ("δ",))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
