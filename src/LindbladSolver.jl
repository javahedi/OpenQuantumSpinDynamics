module LindbladSolver

    using SparseArrays
    using LinearAlgebra
    using DifferentialEquations
    using Expokit

    export LindbladSystem, evolve_Lindblad

struct LindbladSystem
    superop::SparseMatrixCSC{ComplexF64, Int64}
    hamiltonian::SparseMatrixCSC{Float64, Int64}
    lindblad_ops::Vector{SparseMatrixCSC{Float64, Int64}}
end

# Constructor for Lindblad system
function LindbladSystem(hamiltonian::SparseMatrixCSC{Float64, Int64}, 
                        lindblad_ops::Vector{SparseMatrixCSC{Float64, Int64}})
    hbar = 1.0
    
    # Validate dimensions: ensure all Lindblad operators have the same dimension as the Hamiltonian
    dim = size(hamiltonian, 1)
    @assert all(size(op) == (dim, dim) for op in lindblad_ops) "Lindblad operators must have the same dimensions as the Hamiltonian."
    
    # Construct the Lindblad superoperator
    superH = -1im / hbar * (kron(I(dim), hamiltonian) - kron(transpose(hamiltonian), I(dim)))
    
    # The Lindblad superoperator (superL)
    superL = sum([
        kron(conj(op), op) - 0.5 * (kron(I(dim), op' * op) + kron(transpose(op) * conj(op), I(dim))) 
        for op in lindblad_ops
    ])
    
    # Total superoperator
    superop = superH + superL
    
    return LindbladSystem(superop, hamiltonian, lindblad_ops)
end

    function _rho_dot!(dρ, state, superop, time)
        dρ .= superop * state
    end
   
    function _solver_ode(solver::LindbladSystem, ρ0::SparseMatrixCSC{ComplexF64, Int}, 
                        time_points::Vector{Float64}, 
                        observables::Vector{SparseMatrixCSC{Float64, Int}},
                        solver_algorithm=AutoVern8(Vern8()))

        # Validate dimensions
        @assert size(ρ0[:], 1) == size(solver.superop, 1) "Initial state must have compatible dimensions with the superoperator."
        @assert all(size(obs) == size(ρ0) for obs in observables) "Observables must have the same dimensions as the initial state."
     
        
        time_span = (time_points[1],time_points[end])

        #ρ0_vector =  Vector(ρ0[:]) # Vector{Float64}
        # ρ0[:] --> SparseVector{Float64, Int64}
        problem = ODEProblem((dρ, ρ, params, time) -> _rho_dot!(dρ, ρ, solver.superop, time), ρ0[:], time_span)
        solution = solve(problem, solver_algorithm, saveat=time_points, verbose=false)

        outputs = Matrix{Float64}(undef, length(time_points), length(observables))
        for (j, observable) in enumerate(observables)
            for i in 1:length(time_points)
                ρt = reshape(solution.u[i], size(ρ0) )
                outputs[i, j] = real(tr(ρt * observable))
            end
        end
    
        return outputs
    end

    function _propagate_ρt(solver::LindbladSystem, ρ0::SparseMatrixCSC{ComplexF64, Int}, time::Float64)
        #ρ0_vector =  Vector(ρ0[:]) # Vector{Float64}
        # ρ0[:] --> SparseVector{Float64, Int64}
        ρt = expmv(time, solver.superop, Vector(ρ0[:]))
        return reshape(ρt, size(ρ0))
    end

   
    function _solver_expm(solver::LindbladSystem, ρ0::SparseMatrixCSC{ComplexF64, Int}, 
                            time_points::Vector{Float64}, 
                            observables::Vector{SparseMatrixCSC{Float64, Int}})
        # Validate dimensions
        @assert size(ρ0[:], 1) == size(solver.superop, 1) "Initial state must have compatible dimensions with the superoperator."
        @assert all(size(obs) == size(ρ0) for obs in observables) "Observables must have the same dimensions as the initial state."
     
        if time_points[1] !=0.0
            time_points .-=time_points[1]
        end

        outputs = Matrix{Float64}(undef, length(time_points), length(observables))
        for (i, time) in enumerate(time_points)
            ρt = _propagate_ρt(solver, ρ0, time )
            for (j, observable) in enumerate(observables)
                outputs[i, j] = real(tr(ρt * observable))
            end
        end

        return outputs
    end
    
    
    function evolve_Lindblad(solver::LindbladSystem, ρ0::SparseMatrixCSC{Float64, Int}, 
                    time_points::Vector{Float64}, 
                    observables::Vector{SparseMatrixCSC{Float64, Int}}, 
                    method::Symbol=:expm)
                    
        time_span = (time_points[1],time_points[end])
        @assert method in [:expm, :ode] "Method must be either :expm or :ode"
        @assert time_span[1] <= minimum(time_points) && time_span[2] >= maximum(time_points) "time_points must be within time_span."

  

        if method == :expm
            return _solver_expm(solver, complex.(ρ0), time_points, observables)
        elseif method == :ode
            return _solver_ode(solver, complex.(ρ0),  time_points, observables)
        end
    end

end  # End of LiouvillianSolver module