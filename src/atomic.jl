export AtomicMeasure

type AtomicMeasure{T, V}
    v::V # Vector/Tuple of variables
    λ::Vector{T} # The measure is sum λ_i * δ_{support_i}
    support::Vector{Vector{T}} # Elements of the finite support
end
function AtomicMeasure{V, S, T}(vars::V, λ::Vector{S}, support::Vector{Vector{T}})
    AtomicMeasure{promote_type(S, T), V}(vars, λ, support)
end

function Base.show(io::IO, η::AtomicMeasure)
    println(io, "Atomic measure on the variables $(η.v) with $(length(η.λ)) atoms:")
    for (λ, x) in zip(η.λ, η.support)
        println(io, " at $x with weight $λ")
    end
end

function Measure(μ::AtomicMeasure{T}, x::AbstractVector{TT}) where {T, TT}
    Measure{T, monomialtype(TT), monovectype(x)}(μ, x)
end
function Measure{T, MT, MVT}(μ::AtomicMeasure{T}, x) where {T, MT, MVT}
    X = monovec(x)
    sum(μ.λ[i] * dirac(X, μ.v=>μ.support[i]) for i in 1:length(μ.λ))
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
function Base.isapprox(μ::AtomicMeasure, ν::AtomicMeasure; kws...)
    m = length(μ.λ)
    if length(ν.λ) != m
        false
    else
        permcomp((i, j) -> isapprox(μ.λ[i], ν.λ[j]; kws...) && isapprox(μ.support[i], ν.support[j]; kws...), m)
    end
end
