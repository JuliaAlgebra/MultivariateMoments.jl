using Test
using MultivariateMoments

@testset "SymMatrix" begin
    Q = SymMatrix([1, 2, 3], 2)
    symmetric_setindex!(Q, 4, 1, 1)
    @test Q.Q == [4, 2, 3]
    symmetric_setindex!(Q, 5, 1, 2)
    @test Q.Q == [4, 5, 3]
    P = copy(Q)
    @test P.n == 2
    symmetric_setindex!(P, 6, 2, 2)
    @test Q.Q == [4, 5, 3]
    @test P.n == 2
    @test P.Q == [4, 5, 6]
    R = map(i -> i - 1, P)
    @test R.n == 2
    @test R.Q == [3, 4, 5]
end
