@testset "MatMeasure" begin
    Mod.@polyvar x y
    @test_throws ArgumentError matmeasure(Measure([1], [x]), [y])
end

# [HL05] Henrion, D. & Lasserre, J-B.
# Detecting Global Optimality and Extracting Solutions of GloptiPoly 2005

@testset "[HL05] Section 2.3" begin
    Mod.@polyvar x y
    η = AtomicMeasure([x, y], [0.4132, 0.3391, 0.2477], [[1, 2], [2, 2], [2, 3]])
    μ = Measure(η, [x^4, x^3*y, x^2*y^2, x*y^3, y^4, x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1])
    ν = matmeasure(μ, [1, x, y, x^2, x*y, y^2])
    for lrc in (SVDChol(), ShiftChol(1e-14))
        atoms = extractatoms(ν, 1e-4, -1, lrc)
        @test !isnull(atoms)
        @test get(atoms) ≈ η
    end
end

function testelements(X, Y, atol)
    @test length(X) == length(Y)
    for y in Y
        @test any(x -> isapprox(x, y, atol=atol), X)
    end
end
@testset "[HL05] Section 3.3.1" begin
    Mod.@polyvar x y z
    U = [1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0;
         2 0 0 0 0 0; # z^2 = 2z
         0 0 0 1 0 0;
         0 2 0 0 0 0; # y^2 = 2y
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         0 0 2 0 0 0] # x^2 = 2x
    # β will be [z, y, x, y*z, x*z, x*y]
    x = monovec([z, y, x, z^2, y*z, y^2, x*z, x*y, x^2])
    # x*β contains x*y*z, x^2*z, x^2*y which are not present so it show fail
    V = MultivariateMoments.build_system(U', x, sqrt(eps(Float64)))
    @test iszerodimensional(V)
    testelements(V, [[2.0, 2.0, 2.0], [2.0, 2.0, 0.0], [2.0, 0.0, 2.0], [0.0, 2.0, 2.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0], [0.0, 0.0, 0.0]], 1e-11)
end

@testset "[HL05] Section 4" begin
    Mod.@polyvar x y
    s3 = sqrt(1/3)
    η = AtomicMeasure([x, y], [0.25, 0.25, 0.25, 0.25], [[-s3, -s3], [-s3, s3], [s3, -s3], [s3, s3]])
    μ = Measure([1/9,     0,     1/9,     0, 1/9,   0,     0,     0,   0, 1/3,   0, 1/3, 0, 0, 1],
                [x^4, x^3*y, x^2*y^2, x*y^3, y^4, x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1])
    ν = matmeasure(μ, [1, x, y, x^2, x*y, y^2])
    for lrc in (SVDChol(), ShiftChol(1e-16))
        atoms = extractatoms(ν, 1e-16, -1, lrc)
        @test !isnull(atoms)
        @test get(atoms) ≈ η
    end
end

@testset "[LJP17] Example ?" begin
    Mod.@polyvar x[1:2]
    η = AtomicMeasure(x, [2.], [[1., 0.]])
    μ = Measure([2., 0.0, 1e-6],
                monomials(x, 2))
    ν = matmeasure(μ, monomials(x, 1))
    atoms = extractatoms(ν, 1e-5)
    @test !isnull(atoms)
    @test get(atoms) ≈ η
end

@testset "[LJP17] Example ?" begin
    Mod.@polyvar x[1:2]
    η = AtomicMeasure(x, [2.53267], [[1., 0.]])
    ν = matmeasure([2.53267 -0.0 -5.36283e-19; -0.0 -5.36283e-19 -0.0; -5.36283e-19 -0.0 7.44133e-6], monomials(x, 2))
    atoms = extractatoms(ν, 1e-5)
    @test !isnull(atoms)
    @test get(atoms) ≈ η
end

@testset "[LJP17] Example ?" begin
    Mod.@polyvar x[1:2]
    η = AtomicMeasure(x, [14.81065219789812], [[1., 0.22784661679191984]])
    ν = matmeasure([14.8107 3.37426 0.769196 0.175447 0.0400656 0.00917416 0.00211654; 3.37426 0.769196 0.175447 0.0400656 0.00917416 0.00211654 0.000498924; 0.769196 0.175447 0.0400656 0.00917416 0.00211654 0.000498924 0.000125543; 0.175447 0.0400656 0.00917416 0.00211654 0.000498924 0.000125543 3.76638e-5; 0.0400656 0.00917416 0.00211654 0.000498924 0.000125543 3.76638e-5 1.62883e-5; 0.00917416 0.00211654 0.000498924 0.000125543 3.76638e-5 1.62883e-5 1.07817e-5; 0.00211654 0.000498924 0.000125543 3.76638e-5 1.62883e-5 1.07817e-5 1.10658e-5], monomials(x, 6))
    for ranktol in (1e-3, 1e-4, 1e-5)
        atoms = extractatoms(ν, 4e-4)
        @test !isnull(atoms)
        @test get(atoms) ≈ η
    end
end
