module Solvers
    using OpenQuantumSpinDynamics
    using SparseArrays

    export solver_function

    function solver_function(model, Cop, ψ0, observables, TIMEPOINTS, solver_type::Symbol, method::Symbol, n_sample::Int, use_gpu::Bool)
        @assert solver_type in [:KrylovArnoldi, :Lindblad, :swf] "solver_type must be :KrylovArnoldi, :Lindblad, or :swf"

        if solver_type == :KrylovArnoldi
            return _solve_with_KrylovArnoldi(model, ψ0, observables, TIMEPOINTS; method=method)
        elseif solver_type == :Lindblad
            ρ0 = sparse(ψ0 * ψ0')  # Convert initial state to density matrix
            return _solve_with_Lindblad(model, Cop, ρ0, observables, TIMEPOINTS; method=method)
        elseif solver_type == :swf
            return _solve_with_StochasticWavefunction(model, Cop, ψ0, observables, TIMEPOINTS; num_samples=n_sample, use_gpu=use_gpu)
        end
    end

    function _solve_with_KrylovArnoldi(model, ψ0, observables, TIMEPOINTS; method::Symbol=:krylov)
        solver = KrylovArnoldiSystem(model.hamiltonian)
        outputs = evolve_KrylovArnodli(solver, ψ0, TIMEPOINTS, observables, method)
        return outputs
    end

    function _solve_with_Lindblad(model, Cop, ρ0, observables, TIMEPOINTS; method::Symbol=:expm)
        solver = LindbladSystem(model.hamiltonian, Cop)
        outputs = evolve_Lindblad(solver, ρ0, TIMEPOINTS, observables, method)
        return outputs
    end

    function _solve_with_StochasticWavefunction(model, Cop, ψ0, observables, TIMEPOINTS; num_samples::Int=100, use_gpu::Bool=false)
        solver = StochasticWavefunctionSystem(model.hamiltonian, Cop)
        mean_values, std_values = evolve_swf(solver, ψ0, TIMEPOINTS, observables, num_samples; use_gpu=use_gpu)
        return mean_values
    end
end