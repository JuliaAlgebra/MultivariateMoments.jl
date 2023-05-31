function build_system(U::AbstractMatrix, basis::MB.MonomialBasis, ztol, args...)
    # System is
    # y = [U 0] * y
    # where y = x[end:-1:1]
    # which is
    # y = U * β
    m = length(basis)
    r = size(U, 1)
    pivots = [findfirst(j -> U[i, j] != 0, 1:m) for i in 1:r]
    if any(isnothing, pivots)
        keep = map(!isnothing, pivots)
        pivots = pivots[keep]
        r = length(pivots)
        U = U[keep, :]
    end
    monos = basis.monomials
    β = monos[pivots]
    system = [MA.operate(dot, β, U[:, i]) - monos[i] for i in eachindex(monos)]
    filter!(!MP.isconstant, system)
    # Type instability here :(
    if MP.mindegree(monos) == MP.maxdegree(monos) # Homogeneous
        projective_algebraic_set(system, Buchberger(ztol), args...)
    else
        algebraic_set(system, Buchberger(ztol), args...)
    end
end

struct Echelon end

function solve(null::MacaulayNullspace, ::Echelon, args...)
    # If M is multiplied by λ, W is multiplied by √λ
    # so we take √||M|| = √nM
    Z = Matrix(null.matrix')
    rref!(Z, null.accuracy / sqrt(size(Z, 2)))
    #r, vals = solve_system(U', μ.x)
    # TODO determine what is better between rank_check and sqrt(rank_check) here
    return build_system(Z, null.basis, √null.accuracy, args...)
end
