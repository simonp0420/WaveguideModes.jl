using Test
using SafeTestsets

# Use SafeTestsets for isolation of tests
@safetestset "RWG Tests" begin
    include("RWG_test.jl")
end

@safetestset "CWG Tests" begin
    include("CWG_test.jl")
end

@safetestset "Integration Tests" begin
    include("integration_test.jl")
end

