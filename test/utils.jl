function testelements(X, Y, atol)
    @test length(X) == length(Y)
    for y in Y
        @show y
        @test any(x -> isapprox(x, y, atol = atol), X)
    end
end
