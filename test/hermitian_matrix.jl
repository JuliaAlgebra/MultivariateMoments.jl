using Test
using MultivariateMoments

@testset "VectorizedHermitianMatrix" begin
    Q = MultivariateMoments.vectorized_hermitian_matrix(Int, (i, j) -> i == j ? i * 2 - 1 : 2 - im, 2, 1:2)
    @test eltype(Q) == Complex{Int}
    @test 1 == @inferred Q[1, 1]
    @test 2 - im == @inferred Q[1, 2]
    @test 2 + im == @inferred Q[2, 1]
    @test 3 == @inferred Q[2, 2]
    symmetric_setindex!(Q, 4, 1, 1)
    @test Q.Q == [4, 2, 3, -1]
    symmetric_setindex!(Q, 5 - 2im, 1, 2)
    @test Q.Q == [4, 5, 3, -2]
    P = copy(Q)
    @test P.n == 2
    symmetric_setindex!(P, 6, 2, 2)
    @test Q.Q == [4, 5, 3, -2]
    @test P.n == 2
    @test P.Q == [4, 5, 6, -2]
    R = map(i -> i - 1, P)
    @test R.n == 2
    @test R.Q == [3, 4, 5, -2]
    @test P.Q == [4, 5, 6, -2]
    S = map(conj, R)
    @test S.n == 2
    @test S.Q == [3, 4, 5, 2]
    @test_throws ErrorException symmetric_setindex!(S, im, 1, 1)
    M = [1       2 + 3im 4 + 5im
         2 - 3im 6       7 + 8im
         4 - 5im 7 - 8im 9]
    N = MultivariateMoments.vectorized_hermitian_matrix(Int, (i, j) -> M[i, j], 3, 3:-1:1)
    @test Matrix(N) == M[3:-1:1, 3:-1:1]
    symmetric_setindex!(N, 4 - 5im, 3, 1)
    @test Matrix(N) != M[3:-1:1, 3:-1:1]
    M[3, 1] = 4 + 5im
    M[1, 3] = 4 - 5im
    @test Matrix(N) == M[3:-1:1, 3:-1:1]
end
