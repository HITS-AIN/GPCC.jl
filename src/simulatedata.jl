"""
    tobs, yobs, σobs = simulatedata(; σ = 0.75)

Returns toy data in 2 arbitrary bands that are useful for verification and illustrative purposes.
Each of the returned outputs is an array of arrays.


Arguments
================
- `σ` controls Gaussian noise added to simulated fluxes


Returned outputs
================

- `tobs`: Array of arrays of observation times.
- `yobs`: Array of arrays of fluxes.
- `σobs`: Array of error measurements.

An example of how the data are organised is given below.
Note that function [`gpcc`](@ref) expects that the data passed to it
is organised in this exact same way.

# Example

```julia-repl
julia> tobs, yobs, σobs = simulatedata(); # produce synthetic data
julia> typeof(tobs), typeof(yobs), typeof(σobs) # all are arrays of arrays
julia> size(tobs), size(yobs), size(σobs)       # all have length L=2 in this example
julia> length.(tobs) # these lines give us the number of observations per band
julia> length.(yobs)
julia> length.(σobs)
julia> using PyPlot # needs to be indepedently installed.
julia> figure(); title("first band")
julia> errorbar(tobs[1], yobs[1], yerr=σobs[1], fmt="o", label="1st band") # plot data of 1st band
```
"""
function simulatedata()


    rg = MersenneTwister(1)

    #---------------------------------------------------------------------
    # Define GP parameters
    #---------------------------------------------------------------------

    ρ = 3.5 # lengthscale

    delays = [0.0; 2.0]

    α  = [1; 1.5] # scaling coefficients

    b  = [6; 15.0] # offset coefficients


    #---------------------------------------------------------------------
    # Report
    #---------------------------------------------------------------------

    for i in 1:length(delays)

        @printf("\nBand %d\n", i)
        @printf("\t delayed by %.2f\n", delays[i])
        @printf("\t scaled by α[%d]=%.2f\n", i,α[i])
        @printf("\t offset by b[%d]=%.2f\n",i, b[i])
    end

    #---------------------------------------------------------------------
    # Data generation parameters
    #---------------------------------------------------------------------

    N = [60; 50] # number of data items per band

    t = [rand(rg, N[1])*20, [rand(rg, 25)*8; 12.0.+rand(rg, 25)*8]]


    #---------------------------------------------------------------------
    # Define Gaussian process to draw noisy targets
    #---------------------------------------------------------------------

    C = delayedCovariance(matern32, α, delays, ρ, t)

    let

        U, S, V = svd(C)

        C = U * Diagonal(max.(1e-6, abs.(S))) * U'

        makematrixsymmetric!(C)

    end


    #---------------------------------------------------------------------
    # Draw targets and arrange in array
    #---------------------------------------------------------------------

    Y = rand(rg, MvNormal(zeros(sum(N)), C))

    y = Vector{Vector{Float64}}(undef, length(delays))

    mark = 0

    for i in 1:length(delays)

        y[i] = Y[mark+1:mark+N[i]] * α[i] .+ b[i] .+ σ*randn(rg, N[i])

        mark += N[i]

    end

    figure(0) ; cla()

    for i in 1:length(delays)

      plot(t[i], y[i], "o", label = @sprintf("delay = %.3f", delays[i]))

    end

    legend()

    return t, y, [σ*ones(size(yᵢ)) for yᵢ in y]


end
