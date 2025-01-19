module QuantumState

    using SparseArrays

    export AbstractInitialState, NeelState, PolarizedState, construct_state

    abstract type AbstractInitialState end

    struct NeelState <: AbstractInitialState
        L::Int
        direction::String
    end

    struct PolarizedState <: AbstractInitialState
        L::Int
        direction::String
    end

    function construct_state(state::NeelState)
        # Validate direction
        if state.direction ∉ ["z", "x"]
            error("Invalid direction: $(state.direction). Choose 'z' or 'x'.")
        end

        # Define basis states
        up = sparse([1.0, 0.0])
        dn = sparse([0.0, 1.0])

        # Transform to x-basis if needed
        if state.direction == "x"
            up = sqrt(0.5) .* (up + dn)
            dn = sqrt(0.5) .* (up - dn)
        end

        # Construct Néel state
        s = up
        for i in 2:state.L
            if i % 2 == 0
                s = kron(s, dn)
            else
                s = kron(s, up)
            end
        end
        return sparse(s)
    end

    function construct_state(state::PolarizedState)
        # Validate direction
        if state.direction ∉ ["z", "x"]
            error("Invalid direction: $(state.direction). Choose 'z' or 'x'.")
        end

        # Define basis states
        up = sparse([1.0, 0.0])
        dn = sparse([0.0, 1.0])

        # Transform to x-basis if needed
        if state.direction == "x"
            up = sqrt(0.5) .* (up + dn)
            dn = sqrt(0.5) .* (up - dn)
        end

        # Construct polarized state
        s = up
        for i in 2:state.L
            s = kron(s, up)
        end
        return sparse(s)
    end

end  # End of QuantumState module