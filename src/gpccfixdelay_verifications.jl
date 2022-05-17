"""
    minopt, pred, posteriorshiftb = gpcc(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 50, ρmin = 0.1, ρmax = 20.0)

Fit Gaussian Process Cross Correlation (GPCC) model for a given vector of delays.

Data passed to the function are organised as arrays of arrays.
The outer array contains L number of inner arrays where L is the number of bands.
The l-th inner arrays hold the data pertaining to the l-th band.
See [`simulatedata`](@ref) for an example of how data are organised.

Input arguments
===============

- `tarray`: Array of arrays of observation times. There are L number of inner arrays. The l-th array holds the observation times of the l-th band.
- `yarray`: Array of arrays of fluxes. Same structure as `tarray`
- `stdarray`: Array of error measurements. Same structure as `tarray`
- `kernel`: Specifies GP kernel function. Options are GPCC.OU / GPCC.rbf / GPCC.matern32
- `delays`: L-dimensional vector of delays.
- `iterations`: maximum number of iterations done when optimising marginal-likelihood of GP, i.e. optimising hyperparameters.
- `seed`: Random seed controls the random sampling of initial solution.
- `numberofrestarts`: Number of times to repeat optimisation in order to avoid suboptimal solutions due to poor initialisation.
- `initialrandom`: Before optimisation begins, a number of random solutions is sampled and the one with the highest likelihood becomes the starting point for the optimisation.
- `ρmin`: minimum value for lengthscale of GP.
- `ρmax`: maximum value for lengthscale of GP.


Returned arguments
==================
- `minopt`: negative log-likelihood reached when optimising GP hyperparameters.
- `predict`: function for predicting on out-of-sample data.
- `posteriorshiftb`: Gaussian posterior for shift parameters returned as an object of type `MvNormal`.

# Example
```julia-repl
julia> tobs, yobs, σobs = simulatedata(); # produce synthetic data
julia> minopt, pred, posteriorshiftb = gpccfixdelay(tobs, yobs, σobs; kernel = GPCC.rbf, delays = [0.0;2.0;6.0], iterations = 1000) # fit GPCC
julia> trange = collect(-10:0.1:25); # define time interval for predictions
julia> μpred, σpred = pred(trange) # obtain predictions
julia> type(μpred), size(μpred) # predictions are also arrays of arrays, organised just like the data
julia> plot(trange, μpred[1], "b") # plot mean predictions for 1st band
julia> fill_between(trange, μpred[1].+σpred[1], μpred[1].-σpred[1], color="b", alpha=0.3) # plot uncertainties for 1st band
```
"""
function gpcc_verification(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 50, ρmin = 0.1, ρmax = 20.0)

    # Same function as below, but easier name for user to call

    gpccfixdelay(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = seed, numberofrestarts = numberofrestarts, initialrandom = initialrandom, ρmin = ρmin, ρmax = ρmax)


end


function gpccfixdelay(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 50, ρmin = 0.1, ρmax = 20.0)

    #---------------------------------------------------------------------
    # Fix random seed for reproducibility
    #---------------------------------------------------------------------

    rg = MersenneTwister(seed)


    #---------------------------------------------------------------------
    # Set constants
    #---------------------------------------------------------------------

    JITTER = 1e-8

    L = length(tarray)


    #---------------------------------------------------------------------
    # Sort out dimensions
    #---------------------------------------------------------------------

    Narray = length.(tarray)

    @assert(L == length(delays) == length(yarray) == length(tarray) == length(stdarray))


    #---------------------------------------------------------------------
    # Auxiliary matrices
    #---------------------------------------------------------------------

    Y = reduce(vcat, yarray)           # concatenated time series observations

    Sobs = Diagonal(reduce(vcat, stdarray).^2)

    Q  = Qmatrix(Narray)               # matrix for replicating elements

    μb = map(mean, yarray)             # prior mean

    Σb = 100 * diagm(map(var, yarray)) # inflated prior covariance

    B  = Q * Σb * Q'

    b̄  = Q * μb


    #---------------------------------------------------------------------
    function unpack(param)
    #---------------------------------------------------------------------

        @assert(length(param) == L + 1)

        local scale  = makepositive.(param[1+0L:1*L]) .+ 1e-4

        local ρ      = transformbetween(param[L+1], ρmin, ρmax)

        return scale, ρ

    end


    #---------------------------------------------------------------------
    # Covariance matrix verification
    #---------------------------------------------------------------------

    param = [rand(L); rand()]

    # Calculate covarariance matrix like in code
    scale, ρ = unpack(param)
    scale = ones(length(scale))

    K = delayedCovariance(kernel, scale, delays, ρ, tarray)

    KSobsB = K + Sobs + B

    # Calculate covarariance matrix like in paper
    Kalt = Matrix{Matrix{Float64}}(undef, L, L)

    for i in 1:L

        for j in 1:L

            Kalt[i,j] = zeros(length(tarray[i]), length(tarray[j]))

            for n in 1:length(tarray[i])

                for m in 1:length(tarray[j])

                    Kalt[i,j][n,m] = scale[i] * scale[j] * kernel(tarray[i][n]-delays[i], tarray[j][m]-delays[j]; ℓ²=ρ)  + ((i==j) * (n==m) * stdarray[i][n]^2) + (i==j) * Σb[i,i]

                end

            end

        end

    end


    Kalt, K, Sobs, B

end
