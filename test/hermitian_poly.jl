using Test
using MultivariateMoments

# In SumOfSquares, we want the hermitian to be able to hold MOI variables for which
# im * MOI.VariableIndex(i) is not a Complex{MOI.VariableIndex}
# We test this with polynomials as it behaves similarly.

@testset "VectorizedHermitianMatrix with polynomial" begin
    Mod.@polyvar x y
    function _tests(Q)
        for i in 1:2, j in 1:2
            @test (@inferred Q[i, j]) isa eltype(Q)
        end
        @test x == @inferred Q[1, 1]
        @test x + im * y == @inferred Q[1, 2]
        @test x - im * y == @inferred Q[2, 1]
        @test x == @inferred Q[2, 2]
    end
    q = [x, x, x, y]
    Q = VectorizedHermitianMatrix(q, 2)
    @test eltype(Q) == polynomial_type(x * y, Complex{Int})
    _tests(Q)
    R = VectorizedHermitianMatrix{eltype(q), Float64}(q, 2)
    @test eltype(R) == polynomial_type(x * y, Complex{Float64})
    _tests(R)
end
