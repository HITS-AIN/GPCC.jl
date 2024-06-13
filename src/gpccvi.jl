"""
    loglikel, pred, α, postb, ρ = gpcc(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 5, rhomin = 0.1, rhomax = rhomax, verbose = false)

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
- `kernel`: Specifies GP kernel function. Options are GPCC.OU, GPCC.rbf, GPCC.matern32, GPCC.matern52
- `delays`: L-dimensional vector of delays.
- `iterations`: maximum number of iterations done when optimising marginal-likelihood of GP, i.e. optimising hyperparameters.
- `seed`: Random seed controls the random sampling of initial solution.
- `numberofrestarts`: Number of times to repeat optimisation in order to avoid suboptimal solutions due to poor initialisation (default is 1).
- `initialrandom`: Before optimisation begins, a number of random solutions are sampled and the one with the highest likelihood becomes the starting point for the optimisation.
- `rhomin`: minimum value for lengthscale ρ of GP (default 0.1).
- `rhomax`: maximum value for lengthscale ρ of GP.
- `verbose`: true / false (default). If set to `true`, auxiliary messages will be printed out 


Returned arguments
==================
- `loglikel`: log-likelihood reached when optimising GP hyperparameters.
- `predict`: function for predicting on out-of-sample data.
- `α`: coefficients by which the latent Gaussian process is scaled in each band
- `postb`: Gaussian posterior for shift parameters returned as an object of type `MvNormal`.
- `ρ`: length scale of latent Gaussian Process

# Example
```julia-repl
julia> tobs, yobs, σobs, truedelays = simulatedata(); # produce synthetic data
julia> loglikel, pred, α, postb, ρ = gpcc(tobs, yobs, σobs; kernel = GPCC.matern32, delays = truedelays, iterations = 1000);  # fit GPCC
julia> trange = collect(-10:0.1:25); # define time interval for predictions
julia> μpred, σpred = pred(trange) # obtain predictions
julia> type(μpred), size(μpred) # predictions are also arrays of arrays, organised just like the data
julia> plot(trange, μpred[1], "b") # plot mean predictions for 1st band
julia> fill_between(trange, μpred[1].+σpred[1], μpred[1].-σpred[1], color="b", alpha=0.3) # plot uncertainties for 1st band
```
"""
function gpccvi(tarray, yarray, stdarray; kernel = kernel, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 5, rhomin = 0.1, rhomax = 100.0, ρ = NaN, verbose = false)

    # Same function as below, but easier name for user to call

    _gpccvi(tarray, yarray, stdarray; kernel = kernel, iterations = iterations, seed = seed, numberofrestarts = numberofrestarts, initialrandom = initialrandom, ρmin = rhomin, ρmax = rhomax, ρfixed = ρ, verbose = verbose)


end


function _gpccvi(tarray, yarray, stdarray; kernel = kernel, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 5, ρmin = 0.1, ρmax = 20.0, ρfixed = ρfixed, verbose = verbose)

    #---------------------------------------------------------------------
    # Fix random seed for reproducibility
    #---------------------------------------------------------------------

    rg = MersenneTwister(seed)


    #---------------------------------------------------------------------
    # Set constants
    #---------------------------------------------------------------------

    JITTER = 1e-8

    L = length(tarray)

    OPTIMISEρ = isnan(ρfixed) ? :optimise_ρ : :do_not_optimise_ρ


    #---------------------------------------------------------------------
    # Check dimensions
    #---------------------------------------------------------------------

    @assert(L == length(yarray) == length(tarray) == length(stdarray))


    #---------------------------------------------------------------------
    # Auxiliary matrices
    #---------------------------------------------------------------------

    Y = reduce(vcat, yarray)                   # concatenated fluxes

    Q = Qmatrix(length.(tarray))               # matrix for replicating elements

    Sobs = Diagonal(reduce(vcat, stdarray).^2) # observed noise matrix


    μb = map(mean, yarray)             # prior mean

    Σb = 100 * diagm(map(var, yarray)) # inflated prior covariance

    B  = Q * Σb * Q'

    b̄  = Q * μb

    #---------------------------------------------------------------------
    # Let user know what is being run
    #---------------------------------------------------------------------

   

    #---------------------------------------------------------------------
    # Functions for constraining parameters
    #---------------------------------------------------------------------

    # makeα(x) = makepositive(x) + 1e-8

    # makeρ(x) = transformbetween(x, ρmin, ρmax)


    function unpack(param, ::Val{:do_not_optimise_ρ})

        @assert(length(param) == 2L)

        local α =    (param[0L+1:1L])

        local τ = [0;(param[1L+1:2L-1])]

        local ρ = param[2L]

        return α, τ, ρ

    end

    # function unpack(param, ::Val{:optimise_ρ})

    #     @assert(length(param) == L + 1)

    #     local α = makeα.(param[1:1L])

    #     local ρ = makeρ(param[L+1])

    #     return α, ρ

    # end

    unpack(param) = unpack(param, Val(OPTIMISEρ))


    #---------------------------------------------------------------------
    # Define objective as marginal log-likelihood and auxiliaries
    #---------------------------------------------------------------------

    function objective(α, τ, ρ)

        local K = Symmetric(delayedCovariance(kernel, α, τ, ρ, tarray) + Sobs + B)

        return logpdf(MvNormal(b̄, K), Y) 

    end
    
    # helper(p) = objective(unpack(p)...)

    # return VIdiag(helper, 0.01*randn(rg, 2L-1), iterations = iterations,S=300,show_every=1,test_every=50,Stest=1000)

    helper(p) = objective(unpack(p)...)

    fwd(x) = softplus(x)

    bwd(x) = invsoftplus(x)

    f(x) = fwd.(x)#[fwd.(x[1:L]); x[L+1:2L-1]; fwd(x[2L])]
    
    g(x) = bwd.(x)#[bwd.(x[1:L]); x[L+1:2L-1]; bwd(x[2L])]

    elbo = elbofy(2L, 500, helper, transform = f)

    elbohelper(x) = -elbo(x)


    opt = Optim.Options(show_trace = true, iterations = iterations, show_every = 2)

    θ = optimize(elbohelper, [0.01*randn(rg, 2L);0.1*ones(2L)], NelderMead(), opt).minimizer

    MvNormal(θ[1:2L], θ[2L+1:end]), helper, f, g
    
end
