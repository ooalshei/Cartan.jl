using SafeTestsets

@safetestset "builders" begin
    include("builders_tests.jl")
end

@safetestset "involutions tests" begin
    include("involutions_tests.jl")
end

@safetestset "Cartan decomposition tests" begin
    include("decomposition_tests.jl")
end

@safetestset "involutionless Cartan decomposition tests" begin
    include("involutionless_tests.jl")
end

@safetestset "reductive Cartan decomposition tests" begin
    include("reductive_tests.jl")
end

@safetestset "optimizer tests" begin
    include("optimizer_tests.jl")
end

@safetestset "reductive optimizer tests" begin
    include("reductive_optimizer_tests.jl")
end