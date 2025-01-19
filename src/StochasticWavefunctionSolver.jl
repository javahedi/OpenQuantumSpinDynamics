module StochasticWavefunctionSolver

    using LinearAlgebra
    using SparseArrays
    using Statistics

    # Try to load CUDA and CUDA.CUSPARSE, but don't fail if they're not available
    global GPU_AVAILABLE = false
    try
        using CUDA
        using CUDA.CUSPARSE
        global GPU_AVAILABLE = true
        println("CUDA is available. GPU acceleration enabled.")
    catch
        global GPU_AVAILABLE = false
        println("CUDA is not available. Falling back to CPU.")
    end

    export StochasticWavefunctionSystem, evolve_swf

    # Define the solver structure
    struct StochasticWavefunctionSystem
        hamiltonian::SparseMatrixCSC{ComplexF64, Int64}
        lindblad_ops::Vector{SparseMatrixCSC{ComplexF64, Int64}}
        effective_hamiltonian::SparseMatrixCSC{ComplexF64, Int64}
    end

    # Constructor for the solver
    function StochasticWavefunctionSystem(hamiltonian::SparseMatrixCSC{Float64, Int64}, 
                                          lindblad_ops::Vector{SparseMatrixCSC{Float64, Int64}})
        @assert all(size(L) == size(hamiltonian) for L in lindblad_ops) "Lindblad operators must have the same dimensions as the Hamiltonian."
        
        # Convert Hamiltonian and Lindblad operators to ComplexF64
        hamiltonian_complex = SparseMatrixCSC{ComplexF64, Int64}(hamiltonian)
        lindblad_ops_complex = [SparseMatrixCSC{ComplexF64, Int64}(L) for L in lindblad_ops]
        
        # Compute effective Hamiltonian
        effective_hamiltonian = hamiltonian_complex - (1im / 2) * sum(L' * L for L in lindblad_ops_complex)
        
        return StochasticWavefunctionSystem(hamiltonian_complex, lindblad_ops_complex, effective_hamiltonian)
    end

    # Main evolve function with GPU/CPU selection
    function evolve_swf(solver, ψ0, time_points, observables, num_nsample; 
                    use_adaptive=false, tolerance=1e-6, use_gpu=false)
                    
        if use_gpu && GPU_AVAILABLE
            println("GPU acceleration enabled.")
            return evolve_gpu(solver, ψ0, time_points, observables, num_nsample; 
                             use_adaptive=use_adaptive, tolerance=tolerance)
        else
            if use_gpu
                println("GPU acceleration requested but not available. Falling back to CPU.")
            else
                println("Using CPU implementation.")
            end
            return evolve_cpu(solver, ψ0, time_points, observables, num_nsample; 
                             use_adaptive=use_adaptive, tolerance=tolerance)
        end
    end

    # CPU wrapper function
    function evolve_cpu(solver, ψ0, time_points, observables, num_nsample; 
                        use_adaptive=false, tolerance=1e-6)
        # Validate dimensions
        @assert length(ψ0) == size(solver.hamiltonian, 1) "Initial state must have compatible dimensions with the Hamiltonian."
        @assert all(size(obs) == size(solver.hamiltonian) for obs in observables) "Observables must have the same dimensions as the Hamiltonian."

        # Convert ψ0 to Vector{ComplexF64} if it's a SparseVector
        ψ0_dense = Vector{ComplexF64}(ψ0)

        # Convert observables to SparseMatrixCSC{ComplexF64, Int64} if they are Float64
        observables_complex = [SparseMatrixCSC{ComplexF64, Int64}(obs) for obs in observables]

        # Perform stochastic evolution on the CPU
        return _stochastic_evolution(solver, ψ0_dense, time_points, observables_complex, num_nsample)
    end


    # GPU wrapper function (only defined if CUDA is available)
    if GPU_AVAILABLE
        function evolve_gpu(solver, ψ0, time_points, observables, num_nsample; 
                use_adaptive=false, tolerance=1e-6)
            # Convert ψ0 to Vector{ComplexF64} if it's a SparseVector
            ψ0_dense = Vector{ComplexF64}(ψ0)

            # Move data to the GPU and convert to ComplexF64
            ψ0_gpu           = CuArray{ComplexF64}(ψ0_dense)
            hamiltonian_gpu  = CuSparseMatrixCSC{ComplexF64}(solver.hamiltonian)
            lindblad_ops_gpu = [CuSparseMatrixCSC{ComplexF64}(L) for L in solver.lindblad_ops]
            observables_gpu  = [CuSparseMatrixCSC{ComplexF64}(obs) for obs in observables]

            # Compute effective Hamiltonian on the GPU
            effective_hamiltonian_gpu = hamiltonian_gpu - (1im / 2) * sum(L' * L for L in lindblad_ops_gpu)

            # Create a GPU solver
            solver_gpu = StochasticWavefunctionSystem(hamiltonian_gpu, lindblad_ops_gpu, effective_hamiltonian_gpu)

            # Perform stochastic evolution on the GPU
            return _stochastic_evolution(solver_gpu, ψ0_gpu, time_points, observables_gpu, num_nsample)
        end
    end

    # Generic core evolution function
    function _stochastic_evolution_adaptive(solver::StochasticWavefunctionSystem, ψ0::Vector{ComplexF64}, 
                                time_points::Vector{Float64}, 
                                observables::Vector{SparseMatrixCSC{ComplexF64, Int64}}, 
                                num_nsample::Int64; use_adaptive=false, tolerance=1e-6)

        expectation_values = zeros(Float64, length(time_points), length(observables), num_nsample)

        for s in 1:num_nsample
            ψ = copy(ψ0)  # Start with the initial state
            t = time_points[1]
            dt = time_points[2] - time_points[1]

            # Initial expectation values
            for (j, observable) in enumerate(observables)
                expectation_values[1, j, s] = real(dot(ψ, observable * ψ))
            end

            # Time stepping loop
            for i in 2:length(time_points)
                target_time = time_points[i]
                while t < target_time
                    if use_adaptive
                        dt = compute_adaptive_dt(solver, ψ, dt, tolerance)
                    end
                    dt = min(dt, target_time - t)

                    # Compute the probabilities for jumps
                    dps = [real(dt * dot(ψ, L' * L * ψ)) for L in solver.lindblad_ops]
                    dP = sum(dps)

                    # Determine whether a jump occurs
                    u = rand()
                    if dP < u
                        # No jump, evolve under the effective Hamiltonian
                        temp = (I(size(psi0, 1)) - 1im * Heff' * dt) * waves[end]
                        dψ   = (-1im * solver.effective_hamiltonian) * ψ
                        ψ_new = ψ + dt * dψ
                    else
                        # Jump occurs
                        u = rand()
                        Q = cumsum(dps) / dP
                        k = searchsortedfirst(Q, u)
                        ψ_new = solver.lindblad_ops[k] * ψ
                    end

                    # Normalize the wavefunction
                    ψ_new ./= norm(ψ_new)

                    # Update the wavefunction and time
                    ψ = ψ_new
                    t += dt
                end

                # Compute expectation values for the current state
                for (j, observable) in enumerate(observables)
                    expectation_values[i, j, s] = real(dot(ψ, observable * ψ))
                end
            end
        end



        

        # Compute mean and standard deviation across samples
        mean_values = dropdims(mean(expectation_values, dims=3), dims=3)
        std_values  = dropdims(std(expectation_values, dims=3), dims=3)

        return mean_values, std_values
    end

    # Generic adaptive time stepping function
    function compute_adaptive_dt(solver::StochasticWavefunctionSystem, ψ::Vector{ComplexF64}, 
                                dt, tolerance)
        dψ = (-1im * solver.effective_hamiltonian) * ψ
        new_dt = min(dt, tolerance / norm(dψ))
        return new_dt
    end




    function _stochastic_evolution(solver::StochasticWavefunctionSystem, ψ0::Vector{ComplexF64}, 
                                time_points::Vector{Float64}, 
                                observables::Vector{SparseMatrixCSC{ComplexF64, Int64}}, 
                                num_nsample::Int64)

        # Number of time steps
        m  = length(time_points)
        #dt = time_points[2] - time_points[1]

        # Precompute L' * L for each Lindblad operator
        Ldag_L = [L' * L for L in solver.lindblad_ops]

        # Array to store the mean values
        expectation_values = zeros(Float64, m, length(observables), num_nsample)

    
        ψ = Vector{Vector{ComplexF64}}(undef, m)
        
        # Loop over samples
        for s in 1:num_nsample
            ψ[1] = copy(ψ0)

            # Initial expectation values at time 0
            for (j, observable) in enumerate(observables)
                expectation_values[1, j, s] = real(dot(ψ[1]', observable * ψ[1]))
            end


            for i in 2:m
                #Compute the current time step
                dt = time_points[i] - time_points[i-1]
                # Generate a random number in (0,1]
                u = rand()

                # Compute jump probabilities
                dps = [real(dt * (ψ[i-1]' * (Ldag * ψ[i-1]))) for Ldag in Ldag_L]
                dP = sum(dps)

                # Determine whether a jump occurs
                if dP < u
                    # No jump: evolve under the effective Hamiltonian
                    ψ[i] = (I - 1im * solver.effective_hamiltonian' * dt) * ψ[i-1]
                else
                    # Jump occurs: select the jump operator
                    u = rand()
                    Q = cumsum(dps) / dP
                    k = searchsortedfirst(Q, u)
                    ψ[i] = solver.lindblad_ops[k] * ψ[i-1]
                end

                # Normalize the wavefunction and store it
                ψ[i] ./= norm(ψ[i])

                # Compute expectation values for the current state
                for (j, observable) in enumerate(observables)
                    expectation_values[i, j, s] = real(ψ[i]' * observable * ψ[i])
                end
            end

        end


        # Compute mean and standard deviation across samples
        mean_values = dropdims(mean(expectation_values, dims=3), dims=3)
        std_values  = dropdims(std(expectation_values, dims=3), dims=3)


        return mean_values, std_values
    end

end