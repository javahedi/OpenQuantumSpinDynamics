module Coupling

    """
    Example usage: 
    coupling = LongRangeCouplingDisorder(α, L, N)
    println(get_matrix(coupling))
    """
   
    using Combinatorics
    using Distributions
    using StatsBase
    using Base.Iterators


    export AbstractCoupling, LongRangeCouplingDisorder, 
           LongRangeCouplingClean, NearestNeighborCoupling,
           get_matrix, get_N

    # Define the abstract type
    abstract type AbstractCoupling end

    # Define concrete types
    struct LongRangeCouplingDisorder <: AbstractCoupling
        matrix::Matrix{Float64}  # Coupling matrix
        α::Float64               # Power-law exponent
        L::Int                   # Lattice size
        N::Int                   # Number of spins
        reordered::Bool          # Whether indices are reordered
    end

    struct LongRangeCouplingClean <: AbstractCoupling
        matrix::Matrix{Float64}  # Coupling matrix
        α::Float64               # Power-law exponent
        N::Int                   # Number of spins
    end

    struct NearestNeighborCoupling <: AbstractCoupling
        matrix::Matrix{Float64}  # Coupling matrix
        δ::Float64               # Alternating strength
        N::Int                   # Number of spins
    end

    # Constructor for long-range random coupling
    function LongRangeCouplingDisorder(α::Float64, L::Int, N::Int; reordered::Bool=true)
        

        r_index = sample(1:L, N, replace=false) |> sort

        # Reorder indices if requested
        if reordered
            i_index, r_index = reorder_indices(r_index)
        end
        
        # Initialize the coupling matrix Jmn
        Jmn = zeros(Float64, N, N)
        for i in 1:N-1
            for j in i+1:N
                dij = abs(r_index[i] - r_index[j])
                Jmn[i, j] = 1.0 / dij^α             
            end
        end

        Jmn .+= Jmn'
        return LongRangeCouplingDisorder(Jmn, α, L, N, reordered)
    end

    # Constructor for long-range coupling
    function LongRangeCouplingClean(α::Float64, N::Int)
       
        
        # Initialize the coupling matrix Jmn
        Jmn = zeros(Float64, N, N)
        for i in 1:N-1
            for j in i+1:N
                dij = abs(i - j)
                Jmn[i, j] = 1.0 / dij^α
            end
        end

        Jmn .+= Jmn'
        return LongRangeCouplingClean(Jmn, α, N)
    end


    # Constructor for nearest-neighbor coupling
    function NearestNeighborCoupling(δ::Float64, N::Int)
        

        # Initialize the coupling matrix Jmn
        Jmn = zeros(Float64, N, N)
        for n in 1:N-1
            Jmn[n, n+1] = (1.0 + (-1.0)^n * δ)
        end

        Jmn .+= Jmn'
        return NearestNeighborCoupling(Jmn, δ, N)
    end


    # Helper function to reorder indices
    function reorder_indices(r_index, i_index=[])
        if isempty(i_index)
            i_index = 1:length(r_index)
        end

        reordered_r_index = Int[]
        reordered_i_index = Int[]

        while !isempty(r_index)
            min_distance = Inf
            min_indices = nothing

            for i in eachindex(r_index)
                for j in eachindex(r_index)[i+1:end]
                    distance = abs(r_index[i] - r_index[j])
                    if distance < min_distance
                        min_distance = distance
                        min_indices = (i, j)
                    end
                end
            end

            if min_indices !== nothing
                i, j = min_indices
                push!(reordered_r_index, r_index[i], r_index[j])
                push!(reordered_i_index, i_index[i], i_index[j])
                r_index = filter(x -> x ∉ (r_index[i], r_index[j]), r_index)
                i_index = filter(x -> x ∉ (i_index[i], i_index[j]), i_index)
            end
        end

        return (reordered_i_index, reordered_r_index)
    end

    # Common interface functions
    get_matrix(coupling::AbstractCoupling) = coupling.matrix
    get_N(coupling::AbstractCoupling) = coupling.N

end  # End of Coupling module