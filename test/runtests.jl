using Test
using OpenQuantumSpinDynamics


# PauliOps Tests
@testset "PauliOps Tests" begin
    
    N = 3
    paulis = generate_operators(N)  # Now returns a PauliOp struct
    
    # Test the dimensions of the operators
    @test size(paulis.Xop[1]) == (2^N, 2^N)
    @test size(paulis.Zop[2]) == (2^N, 2^N)
    @test size(paulis.Pop[3]) == (2^N, 2^N)
    @test size(paulis.Nop[1]) == (2^N, 2^N)
end

# Coupling Tests
@testset "Coupling Tests" begin

    # Test 1: LongRangeCouplingDisorder
    α = 2.0
    L = 100
    N = 10
    coupling_disorder = LongRangeCouplingDisorder(α, L, N)

    # Test the dimensions of the coupling matrix
    @test size(get_matrix(coupling_disorder)) == (N, N)
    
    # Test 2: LongRangeCouplingClean
    coupling_clean = LongRangeCouplingClean(α, N)

    # Test the dimensions of the coupling matrix
    @test size(get_matrix(coupling_clean)) == (N, N)
  
    # Test if the matrix is symmetric
    @test isapprox(get_matrix(coupling_clean), transpose(get_matrix(coupling_clean)), atol=1e-8)
end



@testset "Setup Tests" begin
    # Test loading parameters from JSON
    params = setup_parameters("configuration_test.json")
    
    @test haskey(params, "lattice_length")
    @test haskey(params, "lattice_size")
    @test haskey(params, "α")
    
    # Initialize system
    hx, hz, Jxy, Jz, TIMEPOINTS, ψ0, Cop, observables = initialize_system(params)
    
    @test params["lattice_length"] == 60
    @test params["lattice_size"] == 6
    @test params["α"] == 1.0
    
end