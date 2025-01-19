
module SpinModels

using SparseArrays
using OpenQuantumSpinDynamics

export SpinModel, model, update_model!

mutable struct SpinModel
    N::Int
    hamiltonian::SparseMatrixCSC{Float64, Int}
    Jxy::Float64
    Jz::Float64
    hx::Vector{Float64}
    hz::Vector{Float64}
    Jmn::AbstractMatrix{Float64}
    bc::Symbol  # Use Symbol for boundary conditions (:pbc or :obc)
end

"""
    model(N::Int, Jxy::Float64, Jz::Float64; kwargs...)

Constructs a `SpinModel` object with the given parameters.

Arguments:
  - `N`: Number of spins.
  - `Jxy`: XY coupling constant.
  - `Jz`: Z coupling constant.

Keyword arguments:
  - `hx`: External magnetic field in the x-direction (default: `zeros(N)`).
  - `hz`: External magnetic field in the z-direction (default: `zeros(N)`).
  - `Jmn`: Coupling matrix (default: `zeros(N, N)`).
  - `pbc`: Periodic boundary condition (default: `false`).

Returns:
  A `SpinModel` instance.
"""
function model(N::Int, Jxy::Float64, Jz::Float64;
                hx::Vector{Float64}=zeros(N), hz::Vector{Float64}=zeros(N),
                Jmn::AbstractMatrix{Float64}=zeros(Float64, N, N), pbc::Bool=false)

    # Validate input dimensions
    if length(hx) != N || length(hz) != N || size(Jmn) != (N, N)
        error("Dimension mismatch: Ensure hx, hz are of size $N and Jmn is $N x $N.")
    end

    # Initialize the Hamiltonian
    hamiltonian = spzeros(Float64, 2^N, 2^N)
    paulis = generate_operators(N)

    # Boundary condition
    bc = pbc ? :pbc : :obc

    # Add local fields
    @views for i in 1:N
        hamiltonian += hx[i] * paulis.Xop[i]
        hamiltonian += hz[i] * paulis.Zop[i]
    end

    # Add interaction terms
    @views for i in 1:N
        for j in i+1:N
            hamiltonian += 2.0 * Jxy * (paulis.Pop[i] * paulis.Nop[j] + paulis.Nop[i] * paulis.Pop[j]) * Jmn[i, j]
            hamiltonian += Jz * (paulis.Zop[i] * paulis.Zop[j]) * Jmn[i, j]
        end
    end

    return SpinModel(N, hamiltonian, Jxy, Jz, hx, hz, Jmn, bc)
end

"""
    update_model!(model::SpinModel, Jxy::Float64, Jz::Float64; kwargs...)

Updates the Hamiltonian of an existing `SpinModel` object with new parameters.

Arguments:
  - `model`: The `SpinModel` instance to update.
  - `Jxy`: XY coupling constant.
  - `Jz`: Z coupling constant.

Keyword arguments:
  - `hx`: New magnetic field in the x-direction.
  - `hz`: New magnetic field in the z-direction.
  - `Jmn`: New coupling matrix.
  - `pbc`: Periodic boundary condition (default: `false`).

Returns:
  Nothing. Modifies the `SpinModel` in place.
"""
function update_model!(model::SpinModel, Jxy::Float64, Jz::Float64;
                        hx::Vector{Float64}=zeros(model.N),  hz::Vector{Float64}=zeros(model.N),
                        Jmn::AbstractMatrix{Float64}=zeros(Float64, model.N, model.N), pbc::Bool=false)

    # Validate input dimensions
    if length(hx) != model.N || length(hz) != model.N || size(Jmn) != (model.N, model.N)
        error("Dimension mismatch: Ensure hx, hz are of size $(model.N) and Jmn is $(model.N) x $(model.N).")
    end

    # Update Hamiltonian
    hamiltonian = spzeros(Float64, 2^model.N, 2^model.N)
    paulis = generate_operators(model.N)

    @views for i in 1:model.N
        hamiltonian += hx[i] * paulis.Xop[i]
        hamiltonian += hz[i] * paulis.Zop[i]
    end

    @views for i in 1:model.N
        for j in i+1:model.N
            hamiltonian += 2.0 * Jxy * (paulis.Pop[i] * paulis.Nop[j] + paulis.Nop[i] * paulis.Pop[j]) * Jmn[i, j]
            hamiltonian +=        Jz * (paulis.Zop[i] * paulis.Zop[j]) * Jmn[i, j]
        end
    end

    model.hamiltonian = hamiltonian
    model.bc = pbc ? :pbc : :obc
end

end  # End of SpinModel module


#=

The @views macro in Julia is used to avoid unnecessary allocations 
when working with array slices. 
It ensures that operations on array slices (e.g., array[1:3]) 
create views instead of copies of the data. 
This can significantly improve performance, especially for large 
arrays or when working with slices repeatedly.


using BenchmarkTools

function without_views(arr)
    sum = 0.0
    for i in 1:length(arr)
        sum += arr[i]
    end
    return sum
end

function with_views(arr)
    sum = 0.0
    @views for i in 1:length(arr)
        sum += arr[i]
    end
    return sum
end

arr = rand(10^6)

@btime without_views($arr)  # Without @views
@btime with_views($arr)     # With @views

=#