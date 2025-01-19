
using OpenQuantumSpinDynamics
using Logging
using BSON
using UUIDs
using SparseArrays
using Distributed

#----------------------------------------------
# Add worker processes if only one process exists
if nprocs() == 1
    addprocs()
end

@everywhere begin
    using OpenQuantumSpinDynamics
end


#----------------------------------------------
function setup_system(config_file::String)
    # Load parameters from configuration file
    params = setup_parameters(config_file)

    # Initialize system
    hx, hz, Jxy, Jz, TIMEPOINTS, ψ0, Cop, observables = initialize_system(params)

    # Create output directory
    output_path = "$(@__DIR__)/DATA_$(params["method"])/alpha$(params["α"])"
    if !isdir(output_path)
        mkpath(output_path)
        @info "Created output directory: $output_path"
    else
        @info "Output directory already exists: $output_path"
    end

    return params, hx, hz, Jxy, Jz, TIMEPOINTS, ψ0, Cop, observables, output_path
end


#----------------------------------------------
# Load parameters and initialize the system
config_path = joinpath(@__DIR__, "configuration.json")
params, hx, hz, Jxy, Jz, TIMEPOINTS, ψ0, Cop, observables, output_path = setup_system(config_path)

#----------------------------------------------
# Run disorder realizations in parallel
results = pmap(d -> disorder(d, hx, hz, Jxy, Jz, TIMEPOINTS, ψ0, Cop, 
                            observables, params, output_path), 1:params["n_disorder"])

#----------------------------------------------
# Summarize results
summarize_results(results)


# inside the main project folder <OpenQuantumSpinDynamics>
# run :
# julia --project=. -t 4 examples/main.jl