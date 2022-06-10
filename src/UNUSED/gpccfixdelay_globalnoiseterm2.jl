"""
    minopt, pred, posterioroffsetvector, scalingcoeff, lengthscale = gpcc(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 5, rhomin = 0.1, rhomax = rhomax)

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
- `kernel`: Specifies GP kernel function. Options are OU(), RBF(), Matern32(), Matern52()
- `delays`: L-dimensional vector of delays.
- `iterations`: maximum number of iterations done when optimising marginal-likelihood of GP, i.e. optimising hyperparameters.
- `seed`: Random seed controls the random sampling of initial solution.
- `numberofrestarts`: Number of times to repeat optimisation in order to avoid suboptimal solutions due to poor initialisation.
- `initialrandom`: Before optimisation begins, a number of random solutions is sampled and the one with the highest likelihood becomes the starting point for the optimisation.
- `rhomin`: minimum value for lengthscale ρ of GP.
- `rhomax`: maximum value for lengthscale ρ of GP.


Returned arguments
==================
- `minopt`: negative log-likelihood reached when optimising GP hyperparameters.
- `predict`: function for predicting on out-of-sample data.
- `posterioroffsetvector`: Gaussian posterior for shift parameters returned as an object of type `MvNormal`.
- `scalingcoeff`: coefficients by which the latent Gaussian process is scaled in each band
- `lengthscale`: length scale of latent Gaussian Process

# Example
```julia-repl
julia> tobs, yobs, σobs = simulatedata(); # produce synthetic data
julia> minopt, pred, posterioroffsetvector, scalingcoeff, lengthscale = gpcc(tobs, yobs, σobs; kernel = RBF(), delays = [0.0;2.0;6.0], iterations = 1000);  # fit GPCC
julia> trange = collect(-10:0.1:25); # define time interval for predictions
julia> μpred, σpred = pred(trange) # obtain predictions
julia> type(μpred), size(μpred) # predictions are also arrays of arrays, organised just like the data
julia> plot(trange, μpred[1], "b") # plot mean predictions for 1st band
julia> fill_between(trange, μpred[1].+σpred[1], μpred[1].-σpred[1], color="b", alpha=0.3) # plot uncertainties for 1st band
```
"""
function gpcc(tarray, yarray, _stdarray_ignore; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 5, rhomin = 0.1, rhomax = rhomax)

    # Same function as below, but easier name for user to call

    gpccfixdelay(tarray, yarray, _stdarray_ignore; kernel = kernel, τ = delays, iterations = iterations, seed = seed, numberofrestarts = numberofrestarts, initialrandom = initialrandom, ρmin = rhomin, ρmax = rhomax)


end


function gpccfixdelay(tarray, yarray, _stdarray_ignore; kernel = kernel, τ = τ, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 5, ρmin = 0.1, ρmax = 20.0)

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
    # Check dimensions
    #---------------------------------------------------------------------

    Narray = length.(tarray)

    @assert(L == length(τ) == length(yarray) == length(tarray))


    #---------------------------------------------------------------------
    # Auxiliary matrices
    #---------------------------------------------------------------------

    Y = reduce(vcat, yarray)           # concatenated time series observations

    Q  = Qmatrix(Narray)               # matrix for replicating elements

    #---------------------------------------------------------------------
    # Let user know what is being run
    #---------------------------------------------------------------------

    informuser(seed = seed, iterations = iterations, numberofrestarts = numberofrestarts,
                initialrandom = initialrandom, JITTER = JITTER, ρmin = ρmin, ρmax = ρmax, Σb = ones(L,L))


    #---------------------------------------------------------------------
    # Functions for constraining parameters
    #---------------------------------------------------------------------

    makeα(x)  = makepositive(x) + 1e-8

    makeρ(x)  = transformbetween(x, ρmin, ρmax)

    makeσ²(x) = makepositive(x)

    function unpack(param)

        @assert(length(param) == 3L + 1)

        local α =   makeα.(param[0L+1:1L])

        local b =         (param[1L+1:2L])

        local σ² = makeσ².(param[2L+1:3L])

        local ρ =    makeρ(param[3L+1])

        return α, b, σ², ρ

    end


    #---------------------------------------------------------------------
    # Define objective as marginal log-likelihood and auxiliaries
    #---------------------------------------------------------------------

    function objective(α, b, σ², ρ)

        local K = delayedCovariance(kernel, α, τ, ρ, tarray)

        local KSobsB = K + Diagonal(Q*σ²)

        makematrixsymmetric!(KSobsB)

        return logpdf(MvNormal(Q*b, KSobsB), Y)

    end


    # joint optimisation

    objective(param) = objective(unpack(param)...)


    # Define negative objective

    negativeobjective(x) = - objective(x)

    # Auxiliary objectives that catch exception and allow optimiser to continue

    # safeobj = safewrapper(objective)

    safenegativeobj = safewrapper(negativeobjective)



    #---------------------------------------------------------------------
    # Define initial values for lengthscale ρ
    #---------------------------------------------------------------------

    initialρvalues = let

        if numberofrestarts == 1 || numberofrestarts == 2

            # pick initial ρ values randomly

            rand(rg, Uniform(ρmin + 1e-3, ρmax - 1e-3), numberofrestarts)

        else

            # initial ρ values defined on grid

            LinRange(ρmin + 1e-3, ρmax - 1e-3, numberofrestarts)

        end

    end


    @printf("\n\tInitial ρ values are:\n")

    map(x -> @printf("\t%f\n", x), initialρvalues)


    #---------------------------------------------------------------------
    # Returns random values for initial scaling vector α
    #---------------------------------------------------------------------

    sampleα() = map(var, yarray) .* (rand(rg, L) * (1.2 - 0.8) .+ 0.8)

    #---------------------------------------------------------------------
    # Returns random values for initial shift vector v
    #---------------------------------------------------------------------

    sampleb() = map(mean, yarray) .* (rand(rg, L) * (1.2 - 0.8) .+ 0.8)

    #---------------------------------------------------------------------
    # Returns random values for initial noise term σ²
    #---------------------------------------------------------------------

    sampleσ²() = 3 * rand(rg, L)

    #---------------------------------------------------------------------
    # Returns random unconstrained solution
    #---------------------------------------------------------------------

    sampleunconstrainedsolution(i) = [invmakepositive.(sampleα());
                                      sampleb();
                                      invmakepositive.(sampleσ²());
                                      invtransformbetween(initialρvalues[i], ρmin, ρmax)]


    #---------------------------------------------------------------------
    # Function below calls optimiser
    #---------------------------------------------------------------------

    function getsolution(i)

        local opt = Optim.Options(show_trace = false, iterations = iterations, show_every = 2, g_tol=1e-6)

        local randomsolutions = [sampleunconstrainedsolution(i) for _ in 1:initialrandom]

        local bestindex = argmin(map(safenegativeobj, randomsolutions))

        local finalresult = optimize(safenegativeobj,randomsolutions[bestindex], NelderMead(), opt)

        return finalresult

    end


    allresults = [getsolution(i) for i in 1:numberofrestarts]

    result     = allresults[argmin([res.minimum for res in allresults])]

    paramopt   = result.minimizer

    @printf("\n\tOverall minimum is %f\n", result.minimum)


    #---------------------------------------------------------------------
    # instantiate learned matrix and observed variance parameter
    #---------------------------------------------------------------------

    @show α, b, σ², ρ = unpack(paramopt)

    K = delayedCovariance(kernel, α, τ, ρ, tarray)

    KSobsB = K + Diagonal(Q*σ²)

    makematrixsymmetric!(KSobsB)




    #---------------------------------------------------------------------
    # Functions for predicting on test data
    #---------------------------------------------------------------------


    function predictTest(ttest::Union{Array{Array{Float64, 1}, 1}, Array{T} where T<:AbstractRange{S} where S<:Real})

        Q✴  = Qmatrix(length.(ttest))


        # dimensions: N × Ntest
        kB✴ = delayedCovariance(kernel, α, τ, ρ, tarray, ttest)

        # Ntest × Ntest
        cB = delayedCovariance(kernel, α, τ, ρ, ttest)

        # full predictive covariance
        Σpred = cB - kB✴' * (KSobsB \ kB✴)

        makematrixsymmetric!(Σpred)

        Σpred = Σpred + JITTER*I

        # predictive mean

        μpred = kB✴' * (KSobsB \ (Y - (Q*b))) + (Q✴ * b)

        return μpred, Σpred

    end



    function predictTest(ttest::Union{AbstractRange{Float64}, Array{Float64,1}})

        Ntest = length(ttest)

        local μpred, Σpred = predictTest([ttest for _ in 1:L])

        # return predictions per "band" and collapse full covariance to diagonal of standard deviations only

        μ_per_band = [μpred[idx] for idx in Iterators.partition(1:L*Ntest, Ntest)]

        σ_per_band = [sqrt.(max.(diag(Σpred)[idx], 1e-6)) for idx in Iterators.partition(1:L*Ntest, Ntest)]

        return μ_per_band, σ_per_band

    end



    function predictTest(ttest::Array{Array{Float64, 1}, 1},
                         ytest::Array{Array{Float64, 1}, 1},
                         σtest::Array{Array{Float64, 1}, 1})

        local μpred, Σpred = predictTest(ttest)

        local Sobs✴ = Diagonal(reduce(vcat, σtest).^2)

        Σpred = Σpred + Sobs✴

        makematrixsymmetric!(Σpred)

        try

            return logpdf(MvNormal(μpred, Σpred), reduce(vcat, ytest))

        catch exception

            if isa(exception, PosDefException)

                local newΣpred = nearestposdef(Σpred; minimumeigenvalue = 1e-6)

                return logpdf(MvNormal(μpred, newΣpred), reduce(vcat, ytest))

            else

                throw(exception)

            end

        end

    end


    # return:
    # • function value returned from optimisation
    # • prediction function
    # • optimised free parameters

    result.minimum, predictTest, (α, b, σ², ρ)
end
