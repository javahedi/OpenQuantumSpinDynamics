
using Plots
using BSON
using Printf
using LaTeXStrings
using PlutoUI
using Measures
using OpenQuantumSpinDynamics

config_path = joinpath(@__DIR__, "configuration.json")

params = setup_parameters(config_path)
hx, hz, Jxy, Jz, timepoints, ψ0, Cop, observables = initialize_system(params)

function magnetization(method::String, α::Float64)
    ave_mx = zeros(length(timepoints))
    ave_mz = zeros(length(timepoints))

    output_folder = "DATA_$method/alpha$α"

    bson_files = filter(endswith(".bson"), readdir(output_folder, join=true))

    for file in bson_files
        BSON.@load file results

        ave_mx .+= results[:,2]
        ave_mz .+= results[:,1]

    end

    ave_mx ./= params["n_disorder"];
    ave_mz ./= params["n_disorder"];

    return ave_mx./params["lattice_size"], ave_mz./params["lattice_size"]
end


method = "arnoldi"
result_α1 = magnetization(method, 1.0);
result_α2 = magnetization(method, 2.0);
result_α3 = magnetization(method, 3.0);


# Plotting
plt = plot(size=(600, 200), dpi=300, xlim=(10^-1,10^2), ylim=(-0.1,1.01),
            bottom_margin=5mm, xscale=:log10,
            top_margin=0mm, left_margin=5mm)  # Initialize plot with custom size


plot!(timepoints, -1.0 .* result_α1[2] , linecolor=:blue, linewidth=2, label=L"α=1.0")
plot!(timepoints, -1.0 .* result_α2[2], linecolor=:red, linewidth=2, label=L"α=2.0")
plot!(timepoints, -1.0 .* result_α3[2], linecolor=:yellow, linewidth=2, label=L"α=3.0")

xlabel!( "Time")
ylabel!(L"\langle m^z_{st} ⟩",fontsize=24)

savefig("$(@__DIR__)/magnetization_$method.pdf")


# inside the main project folder <OpenQuantumSpinDynamics>
# run :
# julia --project=. examples/plot.jl