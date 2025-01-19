module PauliOps

using SparseArrays
using LinearAlgebra

export generate_operators


# Define the Pauli operator structure
struct PauliOp
    Xop::Vector{SparseMatrixCSC{Float64, Int}}
    Zop::Vector{SparseMatrixCSC{Float64, Int}}
    Pop::Vector{SparseMatrixCSC{Float64, Int}}
    Nop::Vector{SparseMatrixCSC{Float64, Int}}
end

# Function to generate Pauli operators for an N-spin system
function generate_operators(N::Int)
    # Define the single-spin Pauli matrices
    σx = sparse([0. 1.; 1. 0.])
    σz = sparse([1. 0.; 0. -1.])
    σp = sparse([0. 1.; 0. 0.])
    σn = sparse([0. 0.; 1. 0.])

    # Create vectors to hold the operators for each spin position
    Xop = Vector{SparseMatrixCSC{Float64, Int}}()
    Zop = Vector{SparseMatrixCSC{Float64, Int}}()
    Pop = Vector{SparseMatrixCSC{Float64, Int}}()
    Nop = Vector{SparseMatrixCSC{Float64, Int}}()

    # Generate tensor products for each operator
    for i in 0:N-1
        leftdim  = 2^i
        rightdim = 2^(N-i-1)
        I_left   = sparse(I, leftdim, leftdim)
        I_right  = sparse(I, rightdim, rightdim)

        push!(Xop, kron(kron(I_left, σx), I_right))
        push!(Zop, kron(kron(I_left, σz), I_right))
        push!(Pop, kron(kron(I_left, σp), I_right))
        push!(Nop, kron(kron(I_left, σn), I_right))
    end

    # Return the PauliOp struct with all the operators
    return PauliOp(Xop, Zop, Pop, Nop)
end

end  # End of PauliOps module
