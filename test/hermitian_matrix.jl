using Test
using MultivariateMoments

@testset "VectorizedHermitianMatrix" begin
    Q = VectorizedHermitianMatrix([1, 2, 3, -1], 2)
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
end
