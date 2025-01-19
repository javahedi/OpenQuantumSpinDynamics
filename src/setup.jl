# src/setup.jl
using OpenQuantumSpinDynamics
using JSON       
using Logging   


function setup_parameters(json_file::String)
    @info "Loading parameters from JSON file: $json_file"
    params = JSON.parsefile(json_file)

    # Define required keys
    required_keys = [
        "lattice_length", "lattice_size", "γ", "J", "Δ","α", "n_disorder", "n_sample",
        "initial_state_type", "initial_state_direction", "magnetization_type",
        "collapse_operator_type", "use_gpu",  "solver_type", "method"
    ]

    # Check if all required keys are present
    for key in required_keys
        if !haskey(params, key)
            error("Missing required parameter in JSON file: '$key'. Please ensure the JSON file contains all required keys.")
        end
    end

    @info "Parameters successfully loaded and validated."
    @info "Solver type: $(params["solver_type"])"
    @info "Solver method: $(params["method"])"
    @info "Number of disorder realizations: $(params["n_disorder"])"
    @info "Number of montecarlo samples: $(params["n_sample"])"
    @info "GPU acceleration: $(params["use_gpu"] ? "enabled" : "disabled")"

    return params
end

function initialize_system(params::Dict{String, Any})
    L = params["lattice_length"]
    N = params["lattice_size"]
    γ = params["γ"]
    J = params["J"]
    Δ = params["Δ"]   
    α = params["α"]
    initial_state_type = params["initial_state_type"] 
    initial_state_direction = params["initial_state_direction"] 
    magnetization_type = params["magnetization_type"]  
    collapse_operator_type = params["collapse_operator_type"]  

    @info "Initializing system with parameters: L=$L, N=$N, γ=$γ, J=$J, α=$α"

    # Initialize magnetic fields
    hx = zeros(Float64, N)
    hz = zeros(Float64, N)
    Jxy, Jz = J, J*Δ

    # Define time points
    TMIN       = 0.1
    TMAX       = 100.0
    TIMEPOINTS = 10 .^ range(log10(TMIN), log10(TMAX); length=100)

    # Construct the initial state based on the specified type and direction
    @info "Constructing initial state: type=$initial_state_type, direction=$initial_state_direction"
    if initial_state_type == "neel"
        ψ0 = construct_state(NeelState(N, initial_state_direction))
    elseif initial_state_type == "polarized"
        ψ0 = construct_state(PolarizedState(N, initial_state_direction))
    else
        error("Invalid initial_state_type: '$initial_state_type'. Choose 'neel' or 'polarized'.")
    end

    # Create Pauli operators for the system
    @info "Creating Pauli operators for system with $N spins."
    paulis = generate_operators(N)

    # Construct collapse operators based on the specified type
    @info "Constructing collapse operators: type=$collapse_operator_type"
    if collapse_operator_type == "decay"
        Cop = [√(0.5γ) * paulis.Nop[i] for i in 1:N]  # decay
    elseif collapse_operator_type == "dephasing"
        Cop = [√(0.5γ) * paulis.Zop[i] for i in 1:N]  # dephasing
    else
        error("Invalid collapse_operator_type: '$collapse_operator_type'. Choose 'decay' or 'dephasing'.")
    end

    # Construct magnetization operators based on the specified type
    @info "Constructing magnetization operators: type=$magnetization_type"
    if magnetization_type == "staggered"
        Z_obser = sum([(-1)^i .* paulis.Zop[i] for i in 1:N])
        X_obser = sum([(-1)^i .* paulis.Xop[i] for i in 1:N])
    elseif magnetization_type == "uniform"
        Z_obser = sum([paulis.Zop[i] for i in 1:N])
        X_obser = sum([paulis.Xop[i] for i in 1:N])
    else
        error("Invalid magnetization_type: '$magnetization_type'. Choose 'staggered' or 'uniform'.")
    end
    observables = [Z_obser, X_obser]
    @info "System initialization complete."
    return hx, hz, Jxy, Jz, TIMEPOINTS, ψ0, Cop, observables
end