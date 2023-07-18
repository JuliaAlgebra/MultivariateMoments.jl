"""
    struct BorderBasis{T,MT<:AbstractMatrix{T},BT}
        matrix::MT
        standard::BT
        border::BT
    end

This matrix with rows indexed by `standard` and columns indexed by `border`
a `standard`-border basis of the ideal `border .- matrix' * standard`
[LLR08, Section 2.5].
The basis `border` is a superset of *corners* of `standard` and a subset
of the *border* of `standard`.

[LLR08] Lasserre, Jean Bernard and Laurent, Monique, and Rostalski, Philipp.
*Semidefinite characterization and computation of zero-dimensional real radical ideals.*
Foundations of Computational Mathematics 8 (2008): 607-647.
"""
struct BorderBasis{T,MT<:AbstractMatrix{T},BT}
    matrix::MT
    standard::BT
    border::BT
end

function solve(
    b::BorderBasis{T},
    solver::SemialgebraicSets.AbstractMultiplicationMatricesSolver =
    MultivariateMoments.SemialgebraicSets.ReorderedSchurMultiplicationMatricesSolver{T}(),
) where {T}
    vars = MP.variables(b.standard)
    m = length(b.standard)
    # If a monomial `border` is not in `b.border` and it is not a corner
    # then it can be divided by another monomial `x^α in b.border`.
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
    completed_border = Dict{eltype(b.border.monomials),Vector{T}}()
    function known_border_coefficients(border)
        return !isnothing(_index(b.border, border)) || haskey(completed_border, border)
    end
    function border_coefficients(border)
        k = _index(b.border, border)
        if isnothing(k)
            if haskey(completed_border, border)
                return completed_border[border]
            else
                return
            end
        else
            return b.matrix[:, k]
        end
    end
    function try_add_to_border(border, o)
        if known_border_coefficients(border)
            return
        end
        for i in findfirst(1:(o-1))
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
        for std in b.standard.monomials
            for (k, shift) in enumerate(vars)
                if k in view(order, 1:(o-1))
                    continue
                end
                try_add_to_border(shift * std, o)
            end
        end
        for (k, shift) in enumerate(vars)
            if k in view(order, 1:(o-1))
                continue
            end
            if all(b.standard.monomials) do std
                return known_border_coefficients(shift * std)
            end
                o += 1
                order[o] = k
                for (col, std) in enumerate(b.standard.monomials)
                    mult[k][:, col] = border_coefficients(shift * std)
                end
            end
        end
    end
    if o < length(vars)
        # Several things could have gone wrong here:
        # 1) We could be missing corners,
        # 2) We have all corners but we could not complete the border because
        #    there is not topological order working
        # In any case, we don't have all multiplication matrices so we abort
        return
    end
    sols = SS.solve(mult, solver)
    return ZeroDimensionalVariety(sols)
end
