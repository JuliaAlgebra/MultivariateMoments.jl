"""
    struct BorderBasis{D,T,MT<:AbstractMatrix{T},B}
        dependence::BasisDependence{D,B}
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
struct BorderBasis{D,T,MT<:AbstractMatrix{T},B}
    dependence::BasisDependence{D,B}
    matrix::MT
    function BorderBasis(
        dependence::BasisDependence{D,B},
        matrix::AbstractMatrix{T},
    ) where {D,T,B}
        return new{D,T,typeof(matrix),B}(dependence, matrix)
    end
end

function independent_basis(b::BorderBasis{StaircaseDependence})
    return standard_basis(b.dependence, trivial = false)
end
function independent_basis(b::BorderBasis{LinearDependence})
    return independent_basis(b.dependence)
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

BorderBasis{D}(b::BorderBasis{D}) where {D} = b

function BorderBasis{LinearDependence}(b::BorderBasis{StaircaseDependence})
    return BorderBasis(
        convert(BasisDependence{Dependence}, b.dependence),
        b.matrix,
    )
end

function BorderBasis{StaircaseDependence}(b::BorderBasis{LinearDependence})
    d = convert(BasisDependence{StaircaseDependence}, b.dependence)
    rows = _indices(
        independent_basis(b.dependence),
        standard_basis(d, trivial = false),
    )
    cols = _indices(dependent_basis(b.dependence), dependent_basis(d))
    return BorderBasis(d, b.matrix[rows, cols])
end

struct StaircaseSolver{
    T,
    R<:RankCheck,
    M<:SemialgebraicSets.AbstractMultiplicationMatricesSolver,
}
    max_partial_iterations::Int
    max_iterations::Int
    rank_check::R
    solver::M
end
function StaircaseSolver{T}(;
    max_partial_iterations::Int = 0,
    max_iterations::Int = -1,
    rank_check::RankCheck = LeadingRelativeRankTol(Base.rtoldefault(T)),
    solver = SS.ReorderedSchurMultiplicationMatricesSolver{T}(),
) where {T}
    return StaircaseSolver{T,typeof(rank_check),typeof(solver)}(
        max_partial_iterations,
        max_iterations,
        rank_check,
        solver,
    )
end

function solve(
    b::BorderBasis{LinearDependence,T},
    solver::StaircaseSolver = StaircaseSolver{T}(),
) where {T}
    return solve(BorderBasis{StaircaseDependence}(b), solver)
end

function solve(
    b::BorderBasis{<:StaircaseDependence,T},
    solver::StaircaseSolver{T} = StaircaseSolver{T}(),
) where {T}
    d = b.dependence
    dependent = dependent_basis(d)
    vars = MP.variables(d)
    standard = standard_basis(d)
    m = length(standard)
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
        return !isnothing(MB.monomial_index(standard, border)) ||
               !isnothing(MB.monomial_index(dependent, border)) ||
               haskey(completed_border, border)
    end
    function border_coefficients(border)
        k = MB.monomial_index(standard, border)
        if !isnothing(k)
            return SparseArrays.sparsevec([k], [one(T)], m)
        end
        k = MB.monomial_index(dependent, border)
        if !isnothing(k)
            v = zeros(T, m)
            row = 0
            for (i, std) in enumerate(standard.monomials)
                j = MB.monomial_index(d.basis, std)
                if !is_trivial(d.dependence[j])
                    row += 1
                    v[i] = b.matrix[row, k]
                end
            end
            return v
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
        for std in standard.monomials
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
            if all(standard.monomials) do std
                return known_border_coefficients(shift * std)
            end
                o += 1
                order[o] = k
                for (col, std) in enumerate(standard.monomials)
                    mult[k][:, col] = border_coefficients(shift * std)
                end
            end
        end
    end
    @assert o <= length(vars)
    partial = o < length(vars)
    if partial
        # Several things could have gone wrong here:
        # 1) We could be missing corners,
        # 2) We have all corners but we could not complete the border because
        #    there is not topological order working
        # We now try to build new relation by comparing partial multiplication matrices
        # We store them in a vector and reshape in a matrix after as it's easy to append to a vector in-place.
        # a matrix after
        if solver.max_partial_iterations == 0
            return
        else
            com_fix = partial_commutation_fix(
                known_border_coefficients,
                border_coefficients,
                T,
                standard,
                vars,
                solver.rank_check,
            )
        end
    else
        if solver.max_iterations == 0
            Uperp = nothing
        else
            Uperp = commutation_fix(mult, solver.solver.ε)
        end
        com_fix = if isnothing(Uperp)
            nothing
        else
            Uperp, standard
        end
    end
    if isnothing(com_fix)
        # The matrices commute, let's simultaneously diagonalize them
        sols = SS.solve(SS.MultiplicationMatrices(mult), solver.solver)
        return ZeroDimensionalVariety(sols)
    else
        Uperp, Ubasis = com_fix
        # The matrices don't commute, let's find the updated staircase and start again
        new_basis, I1, I2 = MB.merge_bases(Ubasis, dependent)
        new_matrix = Matrix{T}(undef, length(new_basis), size(Uperp, 2))
        I_nontrivial_standard = [
            MB.monomial_index(Ubasis, std) for
            std in standard_basis(b.dependence, trivial = false).monomials
        ]
        Uperpstd = Uperp[I_nontrivial_standard, :]
        for i in axes(new_matrix, 1)
            if iszero(I1[i])
                @assert !iszero(I2[i])
                new_matrix[i, :] = b.matrix[:, I2[i]]' * Uperpstd
            else
                @assert iszero(I2[i])
                new_matrix[i, :] = Uperp[I1[i], :]
            end
        end
        null = MacaulayNullspace(new_matrix, new_basis)
        new_solver = StaircaseSolver{T}(;
            max_partial_iterations = solver.max_partial_iterations - partial,
            max_iterations = solver.max_iterations - !partial,
            solver.rank_check,
            solver.solver,
        )
        return solve(null, ShiftNullspace{StaircaseDependence}(new_solver))
    end
end

function commutation_fix(matrices, ε)
    if isempty(first(matrices))
        return
    end
    leading_F = nothing
    leading_ratio = ε
    for i in eachindex(matrices)
        A = matrices[i]
        for j in eachindex(matrices)
            B = matrices[j]
            # We have 0 = A * B - B * A
            # If it is not zero then left leading eigenvector `u` is a new identity
            # with orthogonal space `Uperp`
            # So if for the Macaulay matrix
            # [M_standard M_border]
            # the nullspace was the column space of
            # [I; B]
            # then for the new Macaulay matrix
            # So if for the Macaulay matrix
            # [M_standard M_border]
            # [u          0       ]
            # the nullspace is the column space of
            # [Uperp; B * Uperp]
            F = LinearAlgebra.svd(A * B - B * A, full = true)
            ratio = first(F.S) / (norm(A) * norm(B))
            if ratio > leading_ratio
                leading_F = F
                leading_ratio = ratio
            end
        end
    end
    if isnothing(leading_F)
        return
    else
        return leading_F.U[:, 2:end]
    end
end

function partial_commutation_fix(
    known_border_coefficients,
    border_coefficients,
    ::Type{T},
    standard::MB.SubBasis{MB.Monomial},
    vars,
    rank_check::RankCheck,
) where {T}
    function shifted_border_coefficients(mono, shift)
        coef = border_coefficients(mono)
        ret = zero(coef)
        unknown = zero(MP.polynomial_type(mono, T))
        for i in eachindex(coef)
            if iszero(coef)
                continue
            end
            shifted = shift * standard.monomials[i]
            j = MB.monomial_index(standard, shifted)
            if !isnothing(j)
                ret[j] += coef[i]
            elseif known_border_coefficients(shifted)
                ret .+= coef[i] .* border_coefficients(shifted)
            else
                MA.operate!(MA.add_mul, unknown, coef[i], shifted)
            end
        end
        return ret, unknown
    end
    new_relations = T[]
    unknowns = MP.polynomial_type(prod(vars), T)[]
    for std in standard.monomials
        for x in vars
            mono_x = x * std
            if !known_border_coefficients(mono_x)
                # FIXME what do we do if one of the two monomials is unknown
                # but the other one is known ?
                continue
            end
            for y in vars
                mono_y = y * std
                if !known_border_coefficients(mono_y)
                    # FIXME what do we do if one of the two monomials is unknown
                    # but the other one is known ?
                    continue
                end
                if isnothing(MB.monomial_index(standard, mono_x))
                    if isnothing(MB.monomial_index(standard, mono_y))
                        coef_xy, unknowns_xy =
                            shifted_border_coefficients(mono_y, x)
                    else
                        mono_xy = x * mono_y
                        if known_border_coefficients(mono_xy)
                            coef_xy = border_coefficients(mono_xy)
                            unknowns_yx = zero(PT)
                        else
                            coef_xy = zeros(length(standard.monomials))
                            unknowns_xy = mono_xy
                        end
                    end
                    coef_yx, unknowns_yx =
                        shifted_border_coefficients(mono_x, y)
                else
                    if !isnothing(MB.monomial_index(standard, mono_y))
                        # Let `f` be `known_border_coefficients`.
                        # They are both standard so we'll get
                        # `f(y * mono_x) - f(x * mono_y)`
                        # which will give a zero column, let's just ignore it
                        continue
                    end
                    mono_yx = y * mono_x
                    if known_border_coefficients(mono_yx)
                        coef_yx = border_coefficients(mono_yx)
                        unknowns_yx = zero(PT)
                    else
                        coef_yx = zeros(length(standard.monomials))
                        unknowns_yx = mono_yx
                    end
                    coef_xy, unknowns_xy =
                        shifted_border_coefficients(mono_y, x)
                end
                append!(new_relations, coef_xy - coef_yx)
                push!(unknowns, unknowns_xy - unknowns_yx)
            end
        end
    end
    standard_part = reshape(
        new_relations,
        length(standard.monomials),
        div(length(new_relations), length(standard.monomials)),
    )
    unknown_monos = MP.merge_monomial_vectors(MP.monomials.(unknowns))
    unknown_part = Matrix{T}(undef, length(unknown_monos), length(unknowns))
    for i in eachindex(unknowns)
        unknown_part[:, i] = MP.coefficients(unknowns[i], unknown_monos)
    end
    basis, I1, I2 = MB.merge_bases(standard, MB.SubBasis{MB.Monomial}(unknown_monos))
    M = Matrix{T}(undef, length(basis.monomials), size(standard_part, 2))
    for i in eachindex(basis.monomials)
        if iszero(I1[i])
            @assert !iszero(I2[i])
            M[i, :] = unknown_part[I2[i], :]
        else
            @assert iszero(I2[i])
            M[i, :] = standard_part[I1[i], :]
        end
    end
    F = LinearAlgebra.svd(M, full = true)
    r = rank_from_singular_values(F.S, rank_check)
    return F.U[:, (r+1):end], basis
end

"""
    Base.@kwdef struct AlgebraicBorderSolver{
        D,
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
    D,
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

function solve(b::BorderBasis{E}, solver::AlgebraicBorderSolver{D}) where {D,E}
    if !(E <: D)
        solve(BorderBasis{D}(b), solver)
    end
    # Form the system described in [HL05, (8)], note that `b.matrix`
    # is the transpose to the matrix in [HL05, (8)].
    # [HL05] Henrion, D. & Lasserre, J-B.
    # *Detecting Global Optimality and Extracting Solutions of GloptiPoly*
    # 2005
    ind = independent_basis(b)
    dep = dependent_basis(b.dependence)
    system = [
        dep.monomials[col] - MP.polynomial(b.matrix[:, col], ind) for
        col in eachindex(dep.monomials)
    ]
    filter!(!MP.isconstant, system)
    V = if MP.mindegree(b.dependence) == MP.maxdegree(b.dependence)
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
    struct AlgebraicFallbackBorderSolver{
        M<:Union{Nothing,SS.AbstractMultiplicationMatricesSolver},
        S<:AlgebraicBorderSolver,
    }
        multiplication_matrices_solver::M
        algebraic_solver::S
    end

Solve with `multiplication_matrices_solver` and if it fails, falls back to
solving the algebraic system formed by the border basis with `algebraic_solver`.
"""
struct AlgebraicFallbackBorderSolver{
    M<:Union{Nothing,SS.AbstractMultiplicationMatricesSolver},
    S<:AlgebraicBorderSolver,
}
    multiplication_matrices_solver::M
    algebraic_solver::S
end
function AlgebraicFallbackBorderSolver(;
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
    return AlgebraicFallbackBorderSolver(
        multiplication_matrices_solver,
        algebraic_solver,
    )
end

function solve(b::BorderBasis, solver::AlgebraicFallbackBorderSolver)
    sol = solve(b, _some_args(solver.multiplication_matrices_solver)...)
    if isnothing(sol)
        sol = solve(b, solver.algebraic_solver)
    end
    return sol
end

function solve(
    b::BorderBasis{LinearDependence},
    solver::AlgebraicFallbackBorderSolver{
        <:Union{Nothing,SS.AbstractMultiplicationMatricesSolver},
        <:AlgebraicBorderSolver{StaircaseDependence},
    },
)
    # We will need to convert in both cases anyway
    return solve(BorderBasis{StaircaseDependence}(b), solver)
end
