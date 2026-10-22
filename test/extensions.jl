module TestExtensions

using Test
import DynamicPolynomials as DP
import TypedPolynomials as TP
import MultivariateMoments as MM
import MathOptInterface as MOI

function test_without_indices(x)
    μ = MM.moment_vector([1.0, 2.0], [one(x), x])
    ν = MM.moment_matrix([1.0 2.0; 2.0 4.0], [one(x), x])
    for op in (MOI.Utilities.map_indices, MOI.Utilities.substitute_variables)
        for m in (μ, ν), container in (m, (m, m), [m], [(m, m)], typeof(m)[])
            @test (@inferred op(error, container)) === container
        end
        for container in ((), Tuple{}[], (μ, ν), [(μ, ν)])
            @test (@inferred op(error, container)) === container
        end
        μ_index = MM.moment_vector([MOI.VariableIndex(1)], [one(x)])
        ν_index = MM.moment_matrix(fill(MOI.VariableIndex(1), 1, 1), [one(x)])
        @test !applicable(op, error, μ_index)
        @test !applicable(op, error, ν_index)
    end
    return
end

@testset "MOI extension" begin
    @test Base.get_extension(MM, :MultivariateMomentsMathOptInterfaceExt) !==
          nothing
    DP.@polyvar x
    TP.@polyvar y
    test_without_indices(x)
    test_without_indices(y)
end

end
