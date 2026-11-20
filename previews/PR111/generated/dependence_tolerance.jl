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

const CANDIDATES = [
    "σ[r+1]/σ[1]/√n" => (σ, r, n, _) -> max(σ[r+1] / σ[1], eps()) / sqrt(n),
    "σ[r+1]/σ[r]·√n" => (σ, r, n, _) -> max(σ[r+1] / σ[r], eps()) * sqrt(n),
    "accuracy·√n" => (σ, r, n, accuracy) -> accuracy * sqrt(n),
]

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
    # Check that the package uses the last candidate
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

twisted_cubic_atoms =
    [[t, t^2, t^3] for t in big.([-0.9, -0.4, 0.1, 0.5, 0.5001, 0.95])]
table([
    moments("twisted cubic, d = 3", [x, y, z], 3, twisted_cubic_atoms, δ) for
    δ in [1e-12, 1e-8, 1e-6]
])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
