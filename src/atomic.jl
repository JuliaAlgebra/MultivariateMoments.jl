"""
    struct WeightedDiracMeasure{T}
        center::Vector{T}
        weight::T
    end

Represents a weighted dirac measure with centered at `center` and with weight `weight`.
"""
struct WeightedDiracMeasure{T,AT<:AbstractVector{T}} # Cannot subtype AbstractMeasureLike, since it does not contain the variables corresponding to the coefficients of the center
    center::AT
    weight::T
end
function WeightedDiracMeasure(center::AbstractVector{S}, weight::T) where {S,T}
    U = promote_type(S, T)
    return WeightedDiracMeasure(
        convert(AbstractVector{U}, center),
        convert(U, weight),
    )
end

"""
    struct AtomicMeasure{T, AT, V} <: AbstractMeasureLike{T}
        variables::V                           # Vector/Tuple of variables
        atoms::Vector{WeightedDiracMeasure{T, AT}} # Atoms of the measure
    end

An measure is said to be *atomic* if it is a sum of weighted dirac measures.
For instance, ``\\eta = 2 \\delta_{(1, 0)} + 3 \\delta_{(1/2, 1/2)}`` is an atomic measure since it is a sum of the diracs centered at ``(1, 0)`` and ``(1/2, 1/2)`` and weighted respectively by 2 and 3.
That is, ``\\mathbb{E}_{\\eta}[p] = 2p(1, 0) + 3p(1/2, 1/2)``.

The `AtomicMeasure` struct stores an atomic measure where `variables` contains the variables and `atoms` contains atoms of the measure.
"""
struct AtomicMeasure{T,AT,V} <: AbstractMeasureLike{T}
    variables::V                           # Vector/Tuple of variables
    atoms::Vector{WeightedDiracMeasure{T,AT}} # Atoms of the measure
end

function Base.show(io::IO, η::AtomicMeasure)
    print(
        io,
        "Atomic measure on the variables $(join(η.variables, ", ")) with $(length(η.atoms)) atoms:",
    )
    for atom in η.atoms
        println(io)
        print(io, " at $(atom.center) with weight $(atom.weight)")
    end
end

function moment_vector(η::AtomicMeasure, basis::MB.SubBasis{MB.Monomial})
    return sum(
        atom.weight * dirac(basis.monomials, η.variables => atom.center) for atom in η.atoms
    )
end

function moment_vector(η::AtomicMeasure, monos::AbstractVector{<:MP.AbstractTermLike})
    return moment_vector(η, MB.SubBasis{MB.Monomial}(monos))
end

function expectation(η::AtomicMeasure, p::_APL)
    return sum(δ -> δ.weight * p(η.variables => δ.center), η.atoms)
end
expectation(p::_APL, η::AtomicMeasure) = expectation(η, p)

function compare_modulo_permutation(f, m)
    picked = BitSet()
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
    return true
end

function Base.isapprox(
    η1::WeightedDiracMeasure,
    η2::WeightedDiracMeasure;
    kws...,
)
    return isapprox(η1.weight, η2.weight; kws...) &&
           isapprox(η1.center, η2.center; kws...)
end
function Base.isapprox(η1::AtomicMeasure, η2::AtomicMeasure; kws...)
    m = length(η1.atoms)
    return length(η2.atoms) == m && compare_modulo_permutation(
        (i, j) -> isapprox(η1.atoms[i], η2.atoms[j]; kws...),
        m,
    )
end
