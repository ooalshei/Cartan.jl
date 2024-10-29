using SafeTestsets

@safetestset "Pauli operations tests" begin
    include("pauli_operations_tests.jl")
end

@safetestset "involutions tests" begin
    include("involutions_tests.jl")
end

@safetestset "Hamiltonian generation tests" begin
    include("hamiltonian_tests.jl")
end

@safetestset "Cartan decomposition tests" begin
    include("decomposition_tests.jl")
end

@safetestset "involutionless Cartan decomposition tests" begin
    include("involutionless_tests.jl")
end

@safetestset "iterative Cartan decomposition tests" begin
    include("iterative_tests.jl")
end

@safetestset "optimizer tests" begin
    include("optimizer_tests.jl")
end