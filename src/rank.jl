"""
    abstract type RankCheck end

Method for computing the rank with [`rank_from_singular_values`](@ref)
based on a list of singular values.
"""
abstract type RankCheck end

"""
    rank_from_singular_values(σ, check::RankCheck)

Return the rank of a matrix with singular values `σ` (in decreasing order)
using `check`.
"""
function rank_from_singular_values end

"""
    accuracy(σ, r, check::RankCheck)

Returns a value measuring the accuracy of the rank check `check` returning rank
`r`. This is used by [`Echelon`](@ref) to determine the accuracy to use for the
Gaussian elimination.
"""
function accuracy(σ::Vector{T}, r, check::RankCheck) where {T}
    if iszero(r)
        return zero(T)
    elseif r == lastindex(σ)
        return default_accuracy(σ, check)
    else
        return σ[r+1] / σ[1]
    end
end

"""
    doubt(σ, check::RankCheck)

Returns a value measuring the doubt of the rank check `check`.
Lower values means more doubt so less certainty.
This is used by [`FallbackRank`](@ref) for deciding whether the fallback
should be used.
"""
function doubt end

"""
    struct UserRank <: RankCheck
        pagesize::Int
    end

The user chooses the rank given the singular values in a
`REPL.TerminalMenus.RadioMenu` of page size `pagesize`.

## Example

```julia
julia> rank_from_singular_values([1, 1e-1, 5e-2, 1e-5, 5e-6], UserRank())
Choose the last significant singular value:
   1.0
   0.1
 > 0.05
   1.0e-5
   5.0e-6
3
```
"""
struct UserRank <: RankCheck
    pagesize::Int
end
UserRank() = UserRank(8)
import REPL
function rank_from_singular_values(σ, check::UserRank)
    menu = REPL.TerminalMenus.RadioMenu(string.(σ), pagesize = check.pagesize)
    prompt = "Choose the last significant singular value:"
    return REPL.TerminalMenus.request(prompt, menu)
end
default_accuracy(::Vector{T}, ::UserRank) where {T} = Base.rtoldefault(T)

"""
    struct FixedRank <: RankCheck
        r::Int
    end

The rank is hardcoded to `r`, independently of the singular values.

## Example

```jldoctest
julia> rank_from_singular_values([1, 1e-1, 5e-2, 1e-5, 5e-6], FixedRank(3))
3
```
"""
struct FixedRank <: RankCheck
    r::Int
end
rank_from_singular_values(_, check::FixedRank) = check.r
default_accuracy(::Vector{T}, ::FixedRank) where {T} = Base.rtoldefault(T)

"""
    mutable struct FixedRanks <: RankCheck
        r::Vector{Int}
        current::Int
    end

The `i`th rank is hardcoded to `r[i]`, independently of the singular values.
The field `current` indicates how many ranks have already been asked.
When `current` is `length(r)`, no rank can be asked anymore.

## Example

```jldoctest
julia> check = FixedRanks([2, 3])
FixedRanks([2, 3], 0)

julia> rank_from_singular_values([1, 1e-1, 5e-5, 1e-5, 5e-6], check)
2

julia> rank_from_singular_values([1, 1e-1, 5e-2, 1e-5, 5e-6], check)
3
```
"""
mutable struct FixedRanks <: RankCheck
    r::Vector{Int}
    current::Int
end
FixedRanks(r::Vector{Int}) = FixedRanks(r, 0)
function rank_from_singular_values(_, check::FixedRanks)
    if check.current == length(check.r)
        error(
            "Only `$(length(check.r)) ranks have been given to `FixedRanks` so no more rank can be asked anymore.",
        )
    end
    check.current += 1
    return check.r[check.current]
end
default_accuracy(::Vector{T}, ::FixedRanks) where {T} = Base.rtoldefault(T)

function _findfirst(f, σ)
    return something(findfirst(f, eachindex(σ)), lastindex(σ) + 1) - 1
end

"""
    struct AbsoluteRankTol{T} <: RankCheck
        tol::T
    end

The rank is the number of singular values that are strictly larger than `tol`.

## Example

```jldoctest
julia> rank_from_singular_values([1, 1e-1, 5e-5, 1e-5, 5e-6], AbsoluteRankTol(1e-4))
3
```
"""
struct AbsoluteRankTol{T} <: RankCheck
    tol::T
end
function rank_from_singular_values(σ, check::AbsoluteRankTol)
    return _findfirst(i -> σ[i] ≤ check.tol, σ)
end
default_accuracy(σ::Vector, check::AbsoluteRankTol) = check.tol / last(σ)

"""
    struct LeadingRelativeRankTol{T} <: RankCheck
        tol::T
    end

The rank is the number of singular values that are strictly larger than
`tol * maximum(σ)` where `maximum(σ)` is the largest singular value.

## Example

When the matrix is obtained from a homogeneous problem where the scaling
is irrelevant, `LeadingRelativeRankTol` may be preferable to
[`AbsoluteRankTol`](@ref) as shown below
```jldoctest
julia> rank_from_singular_values(1e6 * [1, 1e-1, 5e-2, 1e-5, 5e-6], AbsoluteRankTol(1e-4))
5

julia> rank_from_singular_values(1e6 * [1, 1e-1, 5e-2, 1e-5, 5e-6], LeadingRelativeRankTol(1e-4))
3
```
"""
struct LeadingRelativeRankTol{T} <: RankCheck
    tol::T
end
function rank_from_singular_values(σ, check::LeadingRelativeRankTol)
    return _findfirst(σ) do i
        return σ[i] ≤ check.tol * first(σ)
    end
end
function default_accuracy(σ::Vector, check::LeadingRelativeRankTol)
    return check.tol * first(σ) / last(σ)
end

"""
    struct DifferentialRankTol{T} <: RankCheck
        tol::T
    end

The rank is the number of singular values before a singular value (not included)
is `tol` times the previous one (included).

## Example

It is sometimes difficult to figure out the tolerance to use in
[`LeadingRelativeRankTol`](@ref). For instance, choosing `1e-3` will consider
`1e-3` in the example below as not part of the rank while `DifferentialRankTol`
would include it because it is close to the previous singular value.
```jldoctest
julia> rank_from_singular_values([1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 1e-6, 5e-7], LeadingRelativeRankTol(1e-3))
5

julia> rank_from_singular_values([1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 1e-6, 5e-7], DifferentialRankTol(1e-2))
6
```
"""
struct DifferentialRankTol{T} <: RankCheck
    tol::T
end
function rank_from_singular_values(σ, check::DifferentialRankTol)
    return _findfirst(σ) do i
        if i == 1
            return false
        else
            return σ[i] ≤ check.tol * σ[i-1]
        end
    end
end
default_accuracy(::Vector, check::DifferentialRankTol) = check.tol

"""
    struct LargestDifferentialRank <: RankCheck
    end

The rank is the number of singular values until the singular value that has
the largest ratio with the next singular value.

## Example

It is sometimes difficult to figure out the tolerance to use in
[`DifferentialRankTol`](@ref). For instance, choosing `1e-2` will consider
`1e-2`, `5e-2` and `1e-3` in the example below as not part of the rank while
`LargestDifferentialRank` would include them because there is a largest gap
later.
```jldoctest
julia> rank_from_singular_values([1, 1e-2, 5e-2, 1e-3, 1e-6, 5e-7], DifferentialRankTol(1e-2))
1

julia> rank_from_singular_values([1, 1e-2, 5e-2, 1e-3, 1e-6, 5e-7], LargestDifferentialRank())
4
```
"""
struct LargestDifferentialRank <: RankCheck end
function rank_from_singular_values(σ, ::LargestDifferentialRank)
    return argmax(i -> σ[i] / σ[i+1], 1:(length(σ)-1))
end
function default_accuracy(::Vector{T}, ::LargestDifferentialRank) where {T}
    return Base.rtoldefault(T)
end
function doubt(σ, check::LargestDifferentialRank)
    if length(σ) ≤ 2
        return zero(eltype(σ))
    else
        r = rank_from_singular_values(σ, check)
        best = σ[r+1] / σ[r]
        second_best = minimum(σ[i+1] / σ[i] for i in 1:(length(σ)-1) if i != r)
        return best / second_best
    end
end

"""
    struct FallbackRank{T,D,F} <: RankCheck
        tol::T
        default::D
        fallback::F
    end

Defaults to checking the rank with `default` and falls back to `fallback`
if the [`doubt`](@ref) is strictly larger than `tol`.
By default, `fallback` is [`UserRank`](@ref).

## Example

The advantage of [`UserRank`](@ref) is that the user get to see if the rank
check is ambiguous and act accordingly. The downside is that it might be
cumbersome if there are many rank checks to do. With `FallbackRank`, the user
only has to sort out the tricky ones.
In the example below, the first example is handled by
[`LargestDifferentialRank`](@ref). For the second one, the user sees that this
is a tricky choice and can manually choose one of the two ranks, then see the
result of the rest of his code using this value of the code and then choose
the other rank and see the impact of this different choice.

```jldoctest
julia> check = FallbackRank(1e-1, LargestDifferentialRank())
FallbackRank{Float64, LargestDifferentialRank, UserRank}(0.1, LargestDifferentialRank(), UserRank(8))

julia> rank_from_singular_values([1, 1e-2, 5e-2, 1e-3, 1e-6, 5e-7], check)
4

julia> rank_from_singular_values([1, 1e-2, 5e-2, 1e-3, 5e-6, 5e-7], check)
Choose the last significant singular value:
 > 1.0
   0.01
   0.05
   0.001
   5.0e-6
   5.0e-7
1

julia> rank_from_singular_values([1, 1e-2, 5e-2, 1e-3, 5e-6, 5e-7], check)
Choose the last significant singular value:
   1.0
   0.01
   0.05
 > 0.001
   5.0e-6
   5.0e-7
4
```
"""
struct FallbackRank{T,D<:RankCheck,F<:RankCheck} <: RankCheck
    tol::T
    default::D
    fallback::F
end

FallbackRank(tol, default::RankCheck) = FallbackRank(tol, default, UserRank())

function rank_from_singular_values(σ, check::FallbackRank)
    if doubt(σ, check.default) > check.tol
        return rank_from_singular_values(σ, check.fallback)
    else
        return rank_from_singular_values(σ, check.default)
    end
end

function default_accuracy(σ::Vector{T}, check::FallbackRank) where {T}
    # FIXME We should use `fallback` if the `fallback` was used instead
    return default_accuracy(σ, check.default)
end

"""
    LowRankLDLTAlgorithm

Method for computing a ``r \\times n`` matrix `U` of a ``n \\times n`` rank ``r`` psd matrix `Q` such that `Q = U'U`.
"""
abstract type LowRankLDLTAlgorithm end

"""
    ShiftCholeskyLDLT <: LowRankLDLTAlgorithm

Shift the matrix by `shift` times the identity matrix before cholesky.
"""
struct ShiftCholeskyLDLT{T} <: LowRankLDLTAlgorithm
    shift::T
end

"""
    SVDLDLT <: LowRankLDLTAlgorithm

Use SVD decomposition.
"""
struct SVDLDLT <: LowRankLDLTAlgorithm end

"""
    MultivariateMoments.low_rank_ldlt(Q::AbstractMatrix, dec::LowRankLDLTAlgorithm, ranktol)

Returns a ``r \\times n`` matrix ``U`` of a ``n \\times n`` rank ``r`` positive semidefinite matrix ``Q`` such that ``Q = U^\\top U``.
The rank of ``Q`` is the number of singular values larger than `ranktol```{} \\cdot \\sigma_1`` where ``\\sigma_1`` is the largest singular value.
"""
function low_rank_ldlt end

function low_rank_ldlt(M::AbstractMatrix, algo::LowRankLDLTAlgorithm, tol::Real)
    return low_rank_ldlt(M, algo, LeadingRelativeRankTol(tol))
end

"""
    struct LowRankLDLT{T,Tr,C<:RankCheck}
        L::Matrix{T}
        singular_values::Vector{Tr}
        rank_check::C
    end

Low-Rank cholesky decomposition `L * Diagonal(singular_values) * L'` of size
`(n, r)` of a `n x n` matrix with singular values
`singular_values[1] > ... > singular_values[n]`.
The rank was chosen given `singular_values` using `rank_check` via
the [`rank_from_singular_values`](@ref) function.
"""
struct LowRankLDLT{T,Tr,C<:RankCheck}
    L::Matrix{T}
    singular_values::Vector{Tr}
    rank_check::C
end

"""
    accuracy(chol::LowRankLDLT)

Return the ratio `rtol = σ[r+1]/σ[1]` where `σ` is the vector of singular
values of the matrix for which `chol` is the rank-`r` Cholesky decomposition.
This is a good relative tolerance to use with the matrix as `σ[r+1]` is the
first singular value that was discarded.
"""
function accuracy(chol::LowRankLDLT)
    return accuracy(chol.singular_values, size(chol.L, 2), chol.rank_check)
end

function low_rank_ldlt(
    M::AbstractMatrix,
    dec::ShiftCholeskyLDLT,
    rank_check::RankCheck,
)
    m = LinearAlgebra.checksquare(M)
    U = cholesky(M + dec.shift * I).U
    σs = map(i -> (U[i, i])^2, 1:m)
    J = sortperm(σs, rev = true)
    σ_sorted = σs[J]
    r = rank_from_singular_values(σ_sorted, rank_check)
    return LowRankLDLT(
        U[J[1:r], :]' * Diagonal(inv.(sqrt.(σ_sorted[1:r]))),
        σ_sorted,
        rank_check,
    )
end

function low_rank_ldlt(M::AbstractMatrix, ::SVDLDLT, rank_check::RankCheck)
    F = svd(M)
    r = rank_from_singular_values(F.S, rank_check)
    return LowRankLDLT(F.U[:, 1:r], F.S, rank_check)
end
