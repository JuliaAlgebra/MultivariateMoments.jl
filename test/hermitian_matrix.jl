using Test
using MultivariateMoments

@testset "VectorizedHermitianMatrix" begin
    Q = MultivariateMoments.vectorized_hermitian_matrix(Int, (i, j) -> i == j ? i * 2 - 1 : 2 - im, 2, 1:2)
    @test Q[1:2, 1:2] isa Matrix{Complex{Int}}
    @test Q[1:2, 1:2] == [1 2 - im; 2 + im 3]
    @test square_getindex(Q, 1:2) isa VectorizedHermitianMatrix{Int, Bool, Complex{Int}}
    @test square_getindex(Q, 1:2).Q == [1, 2, 3, -1]
    @test square_getindex(Q, 2:-1:1) isa VectorizedHermitianMatrix{Int, Bool, Complex{Int}}
    @test square_getindex(Q, 2:-1:1).Q == [3, 2, 1, -1]
    @test square_getindex(Q, 1:1) isa VectorizedHermitianMatrix{Int, Bool, Complex{Int}}
    @test square_getindex(Q, 1:1).Q == [1]
    @test square_getindex(Q, 2:2) isa VectorizedHermitianMatrix{Int, Bool, Complex{Int}}
    @test square_getindex(Q, 2:2).Q == [3]
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
    for QQ in [similar(N, (2, 2)), similar(N, 2, 2), similar(typeof(N), (2, 2)), similar(typeof(N), 2, 2)]
        @test typeof(QQ) == typeof(N)
        symmetric_setindex!(QQ, 2, 1, 1)
        @test QQ[1, 1] == 2
        symmetric_setindex!(QQ, 3 + 4im, 1, 2)
        @test QQ[1, 2] == 3 + 4im
        @test QQ[2, 1] == 3 - 4im
        symmetric_setindex!(QQ, 5, 2, 2)
        @test QQ[2, 2] == 5
    end
end
