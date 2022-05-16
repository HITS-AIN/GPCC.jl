"""
    tobs, yobs, σobs = simulatedata(; σ = 0.1, seed = 1, N = [50; 40; 30], ρ = 1.75)

Returns toy data in 3 arbitrary bands that are useful for verification and illustrative purposes.
Each of the 3 returned outputs is an array of arrays.
Each contains L=3 number of inner arrays.

Arguments
=========
- `σ` controls the Gaussian noise added to the simulated data.
- `seed` controls the random seed for generating the simulated data.
- `N` is a 3-dim vector of integers that specifies the number of observations per band.
- `ρ` lengthscale of latent signal drawn from Gaussian process.

Returned outputs
================

- `tobs`: Array of arrays of observation times.
- `yobs`: Array of arrays of fluxes.
- `σobs`: Array of error measurements.

An example of how the data are organised is given below.
Note that function [`gpccfixdelay`](@ref) expects that the data passed to it
is organised in this exact same way.

# Example

```julia-repl
julia> tobs, yobs, σobs = simulatedata(); # produce synthetic data
julia> typeof(tobs), typeof(yobs), typeof(σobs) # all are arrays of arrays
julia> size(tobs), size(yobs), size(σobs)       # all have length L=3 in this example
julia> length.(tobs) # these three lines give us the number of observations per band
julia> length.(yobs)
julia> length.(σobs)
julia> using PyPlot # needs to be indepedently installed.
julia> errorbar(tobs[1], yobs[1], yerr=σobs[1], fmt="o", label="1st band") # plot data of 1st band
```
"""
function simulatedata(; σ = 0.1, seed = 1, N = [50; 40; 30], ρ = 1.75)


    rg = MersenneTwister(seed)

    #---------------------------------------------------------------------
    # Define GP parameters
    #---------------------------------------------------------------------

    @show delays = [0.0; 2.0; 6.0]

    @show scale  = [1; 2; 0.5]

    @show shift  = [5; 6.0; 9]


    #---------------------------------------------------------------------
    #
    #---------------------------------------------------------------------

    aux  = [rand(rg, N[l]).> 0.5 for l in 1:length(delays)]


    function samplefrominterval(a)

        draw(x) = x==0 ? rand(rg, Uniform(-10.0, 6.0)) : rand(rg, Uniform(8.0, 25.0))

        return [draw(aᵢ) for aᵢ in a]

    end


    t = [samplefrominterval(aux[l]) for l in 1:length(delays)]


    #---------------------------------------------------------------------
    # Define Gaussian process to draw noisy targets
    #---------------------------------------------------------------------

    C = delayedCovariance(matern32, scale, delays, ρ, t)

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

        y[i] = Y[mark+1:mark+N[i]] * scale[i] .+ shift[i] .+ σ*randn(rg, N[i])

        mark += N[i]

    end

    figure(0) ; cla()

    for i in 1:length(delays)

      plot(t[i], y[i], "o", label = @sprintf("delay = %.3f", delays[i]))

    end

    legend()

    return t, y, [σ*ones(size(yᵢ)) for yᵢ in y]


end
