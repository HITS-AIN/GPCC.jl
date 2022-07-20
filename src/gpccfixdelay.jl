"""
    minopt, pred, α, b, ρ = gpccb(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 5, rhomin = 0.1, rhomax = rhomax)

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
- `numberofrestarts`: Number of times to repeat optimisation in order to avoid suboptimal solutions due to poor initialisation.
- `initialrandom`: Before optimisation begins, a number of random solutions is sampled and the one with the highest likelihood becomes the starting point for the optimisation.
- `rhomin`: minimum value for lengthscale ρ of GP.
- `rhomax`: maximum value for lengthscale ρ of GP.


Returned arguments
==================
- `minopt`: negative log-likelihood reached when optimising GP hyperparameters.
- `predict`: function for predicting on out-of-sample data.
- `α`: coefficients by which the latent Gaussian process is scaled in each band
- `b`: Gaussian posterior for shift parameters returned as an object of type `MvNormal`.
- `ρ`: length scale of latent Gaussian Process

# Example
```julia-repl
julia> tobs, yobs, σobs, truedelays = simulatedata(); # produce synthetic data
julia> minopt, pred, α, b, ρ = gpccb(tobs, yobs, σobs; kernel = GPCC.matern32, delays = truedelays, iterations = 1000);  # fit GPCC
julia> trange = collect(-10:0.1:25); # define time interval for predictions
julia> μpred, σpred = pred(trange) # obtain predictions
julia> type(μpred), size(μpred) # predictions are also arrays of arrays, organised just like the data
julia> plot(trange, μpred[1], "b") # plot mean predictions for 1st band
julia> fill_between(trange, μpred[1].+σpred[1], μpred[1].-σpred[1], color="b", alpha=0.3) # plot uncertainties for 1st band
```
"""
function gpccb(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 5, rhomin = 0.1, rhomax = rhomax)

    # Same function as below, but easier name for user to call

    gpccbfixdelay(tarray, yarray, stdarray; kernel = kernel, τ = delays, iterations = iterations, seed = seed, numberofrestarts = numberofrestarts, initialrandom = initialrandom, ρmin = rhomin, ρmax = rhomax)


end


function gpccbfixdelay(tarray, yarray, stdarray; kernel = kernel, τ = τ, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 5, ρmin = 0.1, ρmax = 20.0)

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

    @assert(L == length(τ) == length(yarray) == length(tarray) == length(stdarray))


    #---------------------------------------------------------------------
    # Auxiliary matrices
    #---------------------------------------------------------------------

    Y = reduce(vcat, yarray)                   # concatenated fluxes

    Q = Qmatrix(length.(tarray))               # matrix for replicating elements

    Sobs = Diagonal(reduce(vcat, stdarray).^2) # observed noise matrix


    # Closed form solution for shift vector

    b = (Q'Q)\Q'*Y

    Qb = Q*b

    #---------------------------------------------------------------------
    # Let user know what is being run
    #---------------------------------------------------------------------

    informuser(seed = seed, iterations = iterations, numberofrestarts = numberofrestarts,
                initialrandom = initialrandom, JITTER = JITTER, ρmin = ρmin, ρmax = ρmax, Σb = ones(L,L))


    #---------------------------------------------------------------------
    # Functions for constraining parameters
    #---------------------------------------------------------------------

    makeα(x) = makepositive(x) + 1e-8

    makeρ(x) = transformbetween(x, ρmin, ρmax)

    function unpack(param)

        @assert(length(param) == L + 1)

        local α = makeα.(param[1:1L])

        local ρ = makeρ(param[L+1])

        return α, ρ

    end


    #---------------------------------------------------------------------
    # Define objective as marginal log-likelihood and auxiliaries
    #---------------------------------------------------------------------

    function objective(α, ρ)

        local K = delayedCovariance(kernel, α, τ, ρ, tarray) + Sobs

        makematrixsymmetric!(K)

        return logpdf(MvNormal(Qb, K), Y)

    end

    # convenient call

    objective(param) = objective(unpack(param)...)

    # Define negative objective

    negativeobjective(x) = - objective(x)

    # Auxiliary objective catches exceptions

    safenegativeobj = safewrapper(negativeobjective)


    #---------------------------------------------------------------------
    # Define initial values for lengthscale ρ
    #---------------------------------------------------------------------

    initialρvalues = let

        if numberofrestarts == 1 || numberofrestarts == 2

            # pick initial ρ values randomly

            rand(rg, Uniform(ρmin + 1e-3, ρmax - 1e-3), numberofrestarts)

        else

            # initial ρ values on grid

            collect(MiscUtil.logrange(ρmin + 1e-3, ρmax - 1e-3, numberofrestarts))

        end

    end


    @printf("\n\tInitial ρ values are:\n")

    map(x -> @printf("\t%f\n", x), initialρvalues)


    #---------------------------------------------------------------------
    # Returns random values for initial scaling vector α and shift vector v
    #---------------------------------------------------------------------

    sampleα() = map(var, yarray)  .* (rand(rg, L) * (1.2 - 0.8) .+ 0.8)


    #---------------------------------------------------------------------
    # Returns random unconstrained solution
    #---------------------------------------------------------------------

    sampleunconstrainedsolution(i) = [invmakepositive.(sampleα());
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


    #---------------------------------------------------------------------
    # Restart optimisation multiple times as specified in `numberofrestarts`
    #---------------------------------------------------------------------

    allresults = [getsolution(i) for i in 1:numberofrestarts]

    result     = allresults[argmin([res.minimum for res in allresults])]

    paramopt   = result.minimizer

    @printf("\n\tOverall minimum is %f\n", result.minimum)


    #---------------------------------------------------------------------
    # instantiate learned kernel matrix
    #---------------------------------------------------------------------

    @show α, ρ = unpack(paramopt)

    K = delayedCovariance(kernel, α, τ, ρ, tarray)

    KSobsB = K + Sobs

    makematrixsymmetric!(KSobsB)


    #---------------------------------------------------------------------
    # Functions for predicting on test data
    #---------------------------------------------------------------------

    function predictTest(ttest::Union{Array{Array{Float64, 1}, 1}, Array{T} where T<:AbstractRange{S} where S<:Real})

        # matrix for replicating elements
        Q✴ = Qmatrix(length.(ttest))

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

    result.minimum, predictTest, (α, b, ρ)
end
