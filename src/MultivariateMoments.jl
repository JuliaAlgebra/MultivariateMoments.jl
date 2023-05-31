module MultivariateMoments

using LinearAlgebra

import MutableArithmetics as MA

import MultivariatePolynomials as MP

import MultivariateBases as MB

abstract type AbstractMeasureLike{T} end
abstract type AbstractMomentLike{T} <: AbstractMeasureLike{T} end

abstract type AbstractMoment{T} <: AbstractMomentLike{T} end
abstract type AbstractMeasure{T} <: AbstractMeasureLike{T} end

include("moment.jl")
include("measure.jl")
include("expectation.jl")
include("symmatrix.jl")
include("hermitian_matrix.jl")
include("moment_matrix.jl")
include("atomic.jl")

include("rank.jl")
include("extract.jl")
include("echelon.jl")
include("shift.jl")
include("flat.jl")

include("deprecate.jl")

# Taken from JuMP.jl
# This package exports everything except internal symbols, which are defined as those
# whose name starts with an underscore. Macros whose names start with
# underscores are internal as well. If you don't want all of these symbols
# in your environment, then use `import JuMP` instead of `using JuMP`.

const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__; all = true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS ||
       startswith(sym_string, "_") ||
       startswith(sym_string, "@_")
        continue
    end
    if !(
        Base.isidentifier(sym) ||
        (startswith(sym_string, "@") && Base.isidentifier(sym_string[2:end]))
    )
        continue
    end
    @eval export $sym
end

end # module
