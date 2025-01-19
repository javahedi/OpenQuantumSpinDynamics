#__precompile__(false)
module KrylovArnoldiSolver


    using LinearAlgebra
    using SparseArrays
    using Expokit
    include("Arnoldi.jl")
    

    export KrylovArnoldiSystem, evolve_KrylovArnodli
    
    struct KrylovArnoldiSystem
        hamiltonian::SparseMatrixCSC{Float64, Int64}         # Hamiltonian
    end

    
    # Krylov-based time evolution using Expokit's expmv
    function _solver_Krylov(solver::KrylovArnoldiSystem, ψ0::SparseVector{Float64, Int64}, 
                            observables::Vector{SparseMatrixCSC{Float64, Int}},
                            time_points::Vector{Float64})
                            
                             
        ψ = copy(Vector(ψ0)) #expmv function from the Expokit package does not support sparse vectors (SparseVector{Float64, Int64})
        outputs = Matrix{Float64}(undef, length(time_points), length(observables))
        for (j, observable) in enumerate(observables)
            outputs[1, j] = real(ψ' * observable * ψ)
        end

        told = time_points[1]
        for (i, t) in enumerate(time_points[2:end])
            # Use Expokit's expmv for Krylov-based time evolution
            ψ = expmv(-1.0im * (t - told), solver.hamiltonian, ψ)
            told = t
            for (j, observable) in enumerate(observables)
                outputs[i+1, j] = real(ψ' * observable * ψ)
            end
        end
        return outputs
    end

    # Arnoldi-based time evolution
    function _solver_Arnoldi(solver::KrylovArnoldiSystem, ψ0::SparseVector{Float64, Int64}, 
                            observables::Vector{SparseMatrixCSC{Float64, Int}},
                            time_points::Vector{Float64})
        ψ = copy(ψ0)
        
        outputs = Matrix{Float64}(undef, length(time_points), length(observables))
        for (j, observable) in enumerate(observables)
            outputs[1, j] = real(ψ' * observable * ψ)
        end

        told = time_points[1]
        for (i, t) in enumerate(time_points[2:end])
            # Use Expokit's expmv for Krylov-based time evolution
            ψ, _ = arnoldi(solver.hamiltonian, ψ, t - told)
            told = t
            
            for (j, observable) in enumerate(observables)
                outputs[i+1, j] = real(ψ' * observable * ψ)
            end
        end
        return outputs
    end

    
    function evolve_KrylovArnodli(solver::KrylovArnoldiSystem, ψ0::SparseVector{Float64, Int64}, 
                            time_points::Vector{Float64}, 
                            observables::Vector{SparseMatrixCSC{Float64, Int}}, 
                            method::Symbol=:krylov)

        @assert method in [:arnoldi, :krylov] "Method must be either :arnoldi or :krylov"

        if method == :arnoldi
            return _solver_Arnoldi(solver, ψ0, observables, time_points)
        elseif method == :krylov
            return _solver_Krylov(solver, ψ0, observables, time_points)
        end
    end

end  # End of KrylovArnoldiSolver module