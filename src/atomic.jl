export AtomicMeasure, WeightedDiracMeasure

"""
    struct WeightedDiracMeasure{T}
        center::Vector{T}
        weight::T
    end

Represents a weighted dirac measure with centered at `center` and with weight `weight`.
"""
struct WeightedDiracMeasure{T}
    center::Vector{T}
    weight::T
end
function WeightedDiracMeasure(center::Vector{S}, weight::T) where {S, T}
    U = promote_type(S, T)
    WeightedDiracMeasure(convert(Vector{U}, center), convert(U, weight))
end

"""
    struct AtomicMeasure{T, V}
        variables::V                           # Vector/Tuple of variables
        atoms::Vector{WeightedDiracMeasure{T}} # Atoms of the measure
    end

An measure is said to be *atomic* if it is a sum of weighted dirac measures.
For instance, ``\\eta = 2 \\delta_{(1, 0)} + 3 \\delta_{(1/2, 1/2)}`` is an atomic measure since it is a sum of the diracs centered at ``(1, 0)`` and ``(1/2, 1/2)`` and weighted respectively by 2 and 3.
That is, ``\\mathbb{E}_{\\eta}[p] = 2p(1, 0) + 3p(1/2, 1/2)``.

The `AtomicMeasure` struct stores an atomic measure where `variables` contains the variables and `atoms` contains atoms of the measure.
"""
struct AtomicMeasure{T, V}
    variables::V                           # Vector/Tuple of variables
    atoms::Vector{WeightedDiracMeasure{T}} # Atoms of the measure
end

function Base.show(io::IO, η::AtomicMeasure)
    println(io, "Atomic measure on the variables $(η.variables) with $(length(η.atoms)) atoms:")
    for atom in η.atoms
        println(io, " at $(atom.center) with weight $(atom.weight)")
    end
end

measure(η::AtomicMeasure, x) = Measure(η, x)
function Measure(η::AtomicMeasure{T}, x::AbstractVector{TT}) where {T, TT}
    Measure{T, monomialtype(TT), monovectype(x)}(η, x)
end
function Measure{T, MT, MVT}(η::AtomicMeasure{T}, x) where {T, MT, MVT}
    X = monovec(x)
    sum(atom.weight * dirac(X, η.variables => atom.center) for atom in η.atoms)
end

function permcomp(f, m)
    picked = IntSet()
    for i in 1:m
        k = 0
        for j in 1:m
            if !(j in picked) && f(i, j)
                k = j
                break
            end
        end
        if k == 0
            return false
        end
        push!(picked, k)
    end
    true
end

Base.isapprox(η1::WeightedDiracMeasure, η2::WeightedDiracMeasure; kws...) = isapprox(η1.weight, η2.weight; kws...) && isapprox(η1.center, η2.center; kws...)
function Base.isapprox(η1::AtomicMeasure, η2::AtomicMeasure; kws...)
    m = length(η1.atoms)
    if length(η2.atoms) != m
        false
    else
        permcomp((i, j) -> isapprox(η1.atoms[i], η2.atoms[j]; kws...), m)
    end
end
