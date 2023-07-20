"""
    struct BorderBasis{D<:AbstractDependence,T,MT<:AbstractMatrix{T},BT}
        dependence::D
        matrix::MT
    end

This matrix with rows indexed by `standard` and columns indexed by `border`
a `standard`-border basis of the ideal `border .- matrix' * standard`
[LLR08, Section 2.5].
For solving this with a multiplication matrix solver, it is necessary for the
basis `border` to be a superset of the set of *corners* of `standard` and it is
sufficient for it to be the *border* of `standard`.

[LLR08] Lasserre, Jean Bernard and Laurent, Monique, and Rostalski, Philipp.
*Semidefinite characterization and computation of zero-dimensional real radical ideals.*
Foundations of Computational Mathematics 8 (2008): 607-647.
"""
struct BorderBasis{D<:AbstractDependence,T,MT<:AbstractMatrix{T}}
    dependence::D
    matrix::MT
end

function Base.show(io::IO, b::BorderBasis)
    println(io, "BorderBasis with independent rows and dependent columns in:")
    println(io, b.dependence)
    print(io, "And entries in a ", summary(b.matrix))
    isempty(b.matrix) && return
    println(io, ":")
    Base.print_matrix(io, b.matrix)
    return
end

BorderBasis{D}(b::BorderBasis{<:D}) where {D} = b

function BorderBasis{AnyDependence}(b::BorderBasis{<:StaircaseDependence})
    return BorderBasis(convert(AnyDependence, b.dependence), b.matrix)
end

function BorderBasis{StaircaseDependence}(b::BorderBasis{<:AnyDependence})
    d = StaircaseDependence(_ -> true, b.dependent)
    dependent = convert(AnyDependence, d).dependent
    rows = _indices(b.dependence.independent, d.standard)
    cols = _indices(b.dependence.dependent, dependent)
    return BorderBasis(d, b.matrix[rows, cols])
end

function solve(
    b::BorderBasis{<:StaircaseDependence,T},
    solver::SemialgebraicSets.AbstractMultiplicationMatricesSolver = MultivariateMoments.SemialgebraicSets.ReorderedSchurMultiplicationMatricesSolver{
        T,
    }(),
) where {T}
    d = b.dependence
    dependent = convert(AnyDependence, d).dependent
    vars = MP.variables(d.standard)
    m = length(d.standard)
    # If a monomial `border` is not in `dependent` and it is not a corner
    # then it can be divided by another monomial `x^α in dependent`.
    # So there exists a monomial `x^β` such that `border = x^(α + β)`.
    # In particular, there exists a variable `v` different from `shift`
    # that divides `border` and such that `shift * border / v` is in the border.
    # If the multiplication matrix for `v` happens to have been computed
    # before the multiplication matrix for `shift` then we can obtained
    # the column for `shift * border` as the the multiplication matrix for `v`
    # multiplied by the column for `shift * border / v`.
    # Let `V` be the set of variables such that `shift * border / v` is in the
    # border. We need an order such that at least one of the variables of `V`
    # is before `shift`. So we need some kind of topological sort on the
    # directed forward hypergraph with the variables as nodes and the forward hyperedges
    # `shift => V`. Fortunatly, the Kahn's algorithm generalizes itself quite naturally
    # to hypergraphs.
    order = zeros(Int, length(vars))
    mult = Matrix{T}[zeros(T, m, m) for _ in eachindex(vars)]
    completed_border = Dict{eltype(dependent.monomials),Vector{T}}()
    function known_border_coefficients(border)
        return !isnothing(_index(d.standard, border)) ||
               !isnothing(_index(dependent, border)) ||
               haskey(completed_border, border)
    end
    function border_coefficients(border)
        k = _index(d.standard, border)
        if !isnothing(k)
            return SparseArrays.sparsevec([k], [one(T)], m)
        end
        k = _index(dependent, border)
        if !isnothing(k)
            return b.matrix[:, k]
        end
        if haskey(completed_border, border)
            return completed_border[border]
        else
            return
        end
    end
    function try_add_to_border(border, o)
        if known_border_coefficients(border)
            return
        end
        for i in 1:o
            v = vars[order[i]]
            if MP.divides(v, border)
                other_border = MP.div_multiple(border, v)
                a = border_coefficients(other_border)
                if !isnothing(a)
                    completed_border[border] = mult[order[i]] * a
                    return
                end
            end
        end
    end
    # Variation of Kahn's algorithm
    prev_o = -1
    o = 0
    while o > prev_o
        prev_o = o
        # We iterate over `standard` monomials first`
        # and then over variables so that we encounter
        # the monomials `shift * std` in increasing
        # monomial order so that we know that if `try_add_to_border`
        # fails, it will fail again if we run this for loop again with the same `o`.
        # That allows to limit the number of iteration of the outer loop by `length(vars)`
        for std in d.standard.monomials
            for (k, shift) in enumerate(vars)
                if k in view(order, 1:o)
                    continue
                end
                try_add_to_border(shift * std, o)
            end
        end
        for (k, shift) in enumerate(vars)
            if k in view(order, 1:o)
                continue
            end
            if all(d.standard.monomials) do std
                return known_border_coefficients(shift * std)
            end
                o += 1
                order[o] = k
                for (col, std) in enumerate(d.standard.monomials)
                    mult[k][:, col] = border_coefficients(shift * std)
                end
            end
        end
    end
    @assert o <= length(vars)
    if o < length(vars)
        # Several things could have gone wrong here:
        # 1) We could be missing corners,
        # 2) We have all corners but we could not complete the border because
        #    there is not topological order working
        # In any case, we don't have all multiplication matrices so we abort
        return
    end
    sols = SS.solve(SS.MultiplicationMatrices(mult), solver)
    return ZeroDimensionalVariety(sols)
end

"""
    Base.@kwdef struct AlgebraicBorderSolver{
        D<:AbstractDependence,
        A<:Union{Nothing,SS.AbstractGröbnerBasisAlgorithm},
        S<:Union{Nothing,SS.AbstractAlgebraicSolver},
    }
        algorithm::A = nothing
        solver::S = nothing
    end

Solve a border basis using `algorithm` and `solver by first converting it to a
`BorderBasis{D}`.
"""
struct AlgebraicBorderSolver{
    D<:AbstractDependence,
    A<:Union{Nothing,SS.AbstractGröbnerBasisAlgorithm},
    S<:Union{Nothing,SS.AbstractAlgebraicSolver},
}
    algorithm::A
    solver::S
end
function AlgebraicBorderSolver{D}(;
    algorithm::Union{Nothing,SS.AbstractGröbnerBasisAlgorithm} = nothing,
    solver::Union{Nothing,SS.AbstractAlgebraicSolver} = nothing,
) where {D}
    return AlgebraicBorderSolver{D,typeof(algorithm),typeof(solver)}(
        algorithm,
        solver,
    )
end
#function AlgebraicBorderSolver{D}(solver::SS.AbstractAlgebraicSolver) where {D}
#    return AlgebraicBorderSolver{D}(nothing, solver)
#end
#function AlgebraicBorderSolver{D}(algo::SS.AbstractGröbnerBasisAlgorithm) where {D}
#    return AlgebraicBorderSolver{D}(algo)
#end

_some_args(::Nothing) = tuple()
_some_args(arg) = (arg,)

function _combine_deg(combine, deg, a, b)
    if isempty(a)
        return deg(b)
    elseif isempty(b)
        return deg(b)
    else
        return combine(deg(a), deg(b))
    end
end

function _combine_deg(combine, deg, a::MB.MonomialBasis, b::MB.MonomialBasis)
    return _combine_deg(combine, deg, a.monomials, b.monomials)
end

function solve(b::BorderBasis{E}, solver::AlgebraicBorderSolver{D}) where {D,E}
    if !(E <: D)
        solve(BorderBasis{D}(b), solver)
    end
    # Form the system described in [HL05, (8)], note that `b.matrix`
    # is the transpose to the matrix in [HL05, (8)].
    # [HL05] Henrion, D. & Lasserre, J-B.
    # *Detecting Global Optimality and Extracting Solutions of GloptiPoly*
    # 2005
    d = convert(AnyDependence, b.dependence)
    system = [
        d.dependent.monomials[col] -
        MP.polynomial(b.matrix[:, col], d.independent) for
        col in eachindex(d.dependent.monomials)
    ]
    filter!(!MP.isconstant, system)
    mindeg = _combine_deg(min, MP.mindegree, d.independent, d.dependent)
    maxdeg = _combine_deg(max, MP.maxdegree, d.independent, d.dependent)
    V = if mindeg == maxdeg
        # Homogeneous
        projective_algebraic_set(
            system,
            _some_args(solver.algorithm)...,
            _some_args(solver.solver)...,
        )
    else
        algebraic_set(
            system,
            _some_args(solver.algorithm)...,
            _some_args(solver.solver)...,
        )
    end
    return V
end

"""
    struct BorderWithFallback{M<:Union{Nothing,SS.AbstractMultiplicationMatricesSolver},S<:Union{Nothing,SS.AbstractAlgebraicSolver}}
        multiplication_matrices_solver::M
        algebraic_solver::A
    end

Solve with `multiplication_matrices_solver` and if it fails, falls back to
solving the algebraic system formed by the border basis with `algebraic_solver`.
"""
struct BorderWithFallback{
    M<:Union{Nothing,SS.AbstractMultiplicationMatricesSolver},
    S<:AlgebraicBorderSolver,
}
    multiplication_matrices_solver::M
    algebraic_solver::S
end
function BorderWithFallback(;
    multiplication_matrices_solver::Union{
        Nothing,
        SS.AbstractMultiplicationMatricesSolver,
    } = nothing,
    algorithm::Union{Nothing,SS.AbstractGröbnerBasisAlgorithm} = nothing,
    solver::Union{Nothing,SS.AbstractAlgebraicSolver} = nothing,
    algebraic_solver::AlgebraicBorderSolver = AlgebraicBorderSolver{
        StaircaseDependence,
    }(;
        algorithm,
        solver,
    ),
)
    return BorderWithFallback(multiplication_matrices_solver, algebraic_solver)
end

function solve(b::BorderBasis, solver::BorderWithFallback)
    sol = solve(b, _some_args(solver.multiplication_matrices_solver)...)
    if isnothing(sol)
        sol = solve(b, solver.algebraic_solver)
    end
    return sol
end

function solve(
    b::BorderBasis{AnyDependence},
    solver::BorderWithFallback{M,<:AlgebraicBorderSolver{StaircaseDependence}},
) where {M}
    # We will need to convert in both cases anyway
    return solve(BorderBasis{StaircaseDependence}(b), solver)
end
