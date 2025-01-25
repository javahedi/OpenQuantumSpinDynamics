

using Plots
using BSON
using Printf
using LaTeXStrings
using PlutoUI
using Measures
using OpenQuantumSpinDynamics
using Statistics  # For mean and std functions

config_path = joinpath(@__DIR__, "configuration.json")

params = setup_parameters(config_path)
hx, hz, Jxy, Jz, timepoints, ψ0, Cop, observables = initialize_system(params)

function magnetization(method::String, α::Float64)
    # Initialize arrays to store all results for mean and std calculation
    all_mx = Vector{Vector{Float64}}()
    all_mz = Vector{Vector{Float64}}()
    all_xx = Vector{Vector{Float64}}()
    all_zz = Vector{Vector{Float64}}()

    output_folder = "$(@__DIR__)/DATA_$method/alpha$α/"

    bson_files = filter(endswith(".bson"), readdir(output_folder, join=true))

    for file in bson_files
        BSON.@load file results

        # Append results to the respective arrays
        push!(all_mx, results[:,2])
        push!(all_mz, results[:,1])
        push!(all_xx, results[:,4])
        push!(all_zz, results[:,3])
    end

    # Convert to matrices for easier computation
    all_mx = hcat(all_mx...)
    all_mz = hcat(all_mz...)
    all_xx = hcat(all_xx...)
    all_zz = hcat(all_zz...)

    # Calculate mean and std using Statistics library
    ave_mx = mean(all_mx, dims=2)
    ave_mz = mean(all_mz, dims=2)
    ave_xx = mean(all_xx, dims=2)
    ave_zz = mean(all_zz, dims=2)

    std_mx = std(all_mx, dims=2)
    std_mz = std(all_mz, dims=2)
    std_xx = std(all_xx, dims=2)
    std_zz = std(all_zz, dims=2)

    # Normalize by system size
    N = params["lattice_size"]
    ave_mx ./= N
    ave_mz ./= N
    ave_xx ./= (N * (N - 1) / 2.0)
    ave_zz ./= (N * (N - 1) / 2.0)

    std_mx ./= N
    std_mz ./= N
    std_xx ./= (N * (N - 1) / 2.0)
    std_zz ./= (N * (N - 1) / 2.0)

    return ave_mx, ave_mz, ave_xx, ave_zz, std_mx, std_mz, std_xx, std_zz
end

method = params["method"]
#result_α1 = magnetization(method, 1.0);
result_α2 = magnetization(method, 2.0);
#result_α3 = magnetization(method, 3.0);


# Plotting
plt = plot(size=(600, 200), dpi=300, xlim=(10^-1,10^2), ylim=(-0.1,1.01),
            bottom_margin=5mm, xscale=:log10, linewidth=2,
            top_margin=0mm, left_margin=5mm)  # Initialize plot with custom size


#plot!(timepoints, -1.0 .* result_α1[2] , linecolor=:blue, linewidth=2, label=L"α=1.0")
#plot!(timepoints, -1.0 .* result_α2[2], linecolor=:red, linewidth=2, label=L"α=2.0")
#plot!(timepoints, -1.0 .* result_α3[2], linecolor=:blue, linewidth=2, label=L"α=2.0")
#xlabel!( "Time")
#ylabel!(L"\langle m^z_{st} ⟩",fontsize=24)

plot!(timepoints,  -1.0 .* result_α2[2], ribbon=result_α2[6], label=L"α=2.0,~\langle m^z_{st}\rangle", color=:blue)
plot!(timepoints, result_α2[4],  ribbon=result_α2[8], label=L"α=2.0,~\langle s_i^zs_j^z\rangle", color=:red)

xlabel!( "Time")
savefig("$(@__DIR__)/magnetization_$method.pdf")


# inside the main project folder <OpenQuantumSpinDynamics>
# run :
# julia --project=. examples/plot.jl