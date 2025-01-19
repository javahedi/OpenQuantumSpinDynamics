module Disorder
    using OpenQuantumSpinDynamics
    using UUIDs
    using BSON
    using Logging
    using Distributed

    export disorder

    function disorder(d, hx, hz, Jxy, Jz, TIMEPOINTS, ψ0, Cop, observables, params, output_path)
        worker_id = myid()
        try
            # Create coupling for each realization
            coupling = LongRangeCouplingDisorder(params["α"], params["lattice_length"], params["lattice_size"]; reordered=true)

            # Initialize the model
            spin_model = model(params["lattice_size"], Jxy, Jz, hx=hx, hz=hz, Jmn=coupling.matrix, pbc=false)

            # Solve for expectation values
            results = solver_function(
                            spin_model, Cop, ψ0, observables, TIMEPOINTS,
                            Symbol(params["solver_type"]), Symbol(params["method"]),
                            Int(params["n_sample"]), Bool(params["use_gpu"])
                        )

            # Save results with a unique ID
            unique_id = UUIDs.uuid1()
            output_file = joinpath(output_path, "disorder_$(unique_id).bson")
            BSON.@save output_file results

            @info "Worker $worker_id: Disorder realization $d completed successfully."
            return true
        catch e
            @error "Worker $worker_id: Disorder realization $d failed with error: $e" exception = (e, catch_backtrace())
            return false
        end
    end
end