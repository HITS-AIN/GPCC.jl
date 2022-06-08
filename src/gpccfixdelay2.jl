"""
    minopt, pred, posterioroffsetvector, scalingcoeff, lengthscale = gpcc(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 5, rhomin = 0.1, rhomax = 20.0)

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
function gpcc(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 1, rhomin = 0.1, rhomax = 20.0)

    # Same function as below, but easier name for user to call

    gpccfixdelay(tarray, yarray, stdarray; kernel = kernel, τ = delays, iterations = iterations, seed = seed, numberofrestarts = numberofrestarts, initialrandom = initialrandom, ρmin = rhomin, ρmax = rhomax)


end


function gpccfixdelay(tarray, yarray, stdarray; kernel = kernel, τ = τ, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 1, ρmin = 0.1, ρmax = 20.0)

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

    @assert(L == length(τ) == length(yarray) == length(tarray) == length(stdarray))


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
    # Let user know what is being run
    #---------------------------------------------------------------------

    informuser(seed = seed, iterations = iterations, numberofrestarts = numberofrestarts,
                initialrandom = initialrandom, JITTER = JITTER, ρmin = ρmin, ρmax = ρmax, Σb = Σb)


    #---------------------------------------------------------------------
    # Functions for constraining parameters
    #---------------------------------------------------------------------

    makeρ(x) = transformbetween(x, ρmin, ρmax)

    makeα(x) = makepositive(x)

    function unpack(param)

        @assert(length(param) == L + 1)

        local α = makeα.(param[1:L])

        local ρ = makeρ(param[L+1])

        α, ρ

    end


    #---------------------------------------------------------------------
    # Define objective as marginal log-likelihood and auxiliaries
    #---------------------------------------------------------------------

    function objective(α, ρ)

        local K = delayedCovariance(kernel, α, τ, ρ, tarray)

        local KSobsB = K + Sobs + B

        makematrixsymmetric!(KSobsB)

        return logpdf(MvNormal(b̄, KSobsB), Y)

    end

    # optimisation of scalings only

    objective_fixed_lengthscale(paramα, paramρ) = objective(unpack([paramα; paramρ])...)

    # joint optimisation of scalings and lengthscale

    objective(param) = objective(unpack(param)...)


    # Define negative objective

    negativeobjective(x) = - objective(x)

    # Auxiliary objectives that catch exception and allow optimiser to continue

    safeobj = safewrapper(objective)

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
    # Returns random values for initial scalings α
    #---------------------------------------------------------------------

    sampleα() = map(var, yarray) .* (rand(rg, L) * (1.2 - 0.8) .+ 0.8)



    #---------------------------------------------------------------------
    # Function below calls optimiser.
    # Optimisation starts with randomly sampled values for the scalings α
    # and a value on a grid for the lengthscale ρ.
    # Optimisation initially fixed ρ and optimises only ρ.
    # After that joint optimisation follows.
    #---------------------------------------------------------------------

    function getsolution(i)

        # create vectors of initial solutions

        local randomsolutions = [[invmakepositive.(sampleα()); invtransformbetween(initialρvalues[i], ρmin, ρmax)] for _ in 1:initialrandom]

        # evaluate all initial solutions and pick best one

        local bestindex = argmin(map(safeobj, randomsolutions))

        # options for optimiser called below

        local opt = Optim.Options(show_trace = false, iterations = iterations, show_every = 2, g_tol=1e-6)


        # define objective for optimising scalings α only

        local initialparamρ = randomsolutions[bestindex][L+1]

        local initialparamα = randomsolutions[bestindex][1:L]

        local negativeobjective_fixed_lengthscale(x) = - objective_fixed_lengthscale(x, initialparamρ)

        local safenegativeobjective_fixed_lengthscale = safewrapper(negativeobjective_fixed_lengthscale)


        # optimising scalings α only

        @printf("\n\tRestart #%d: starting optimisation with ρ=%f and negative objective %f\n", i, makeρ(initialparamρ), safenegativeobjective_fixed_lengthscale(initialparamα))

        local intermediateresult = optimize(safenegativeobjective_fixed_lengthscale, initialparamα, NelderMead(), opt)

        local optparamα = intermediateresult.minimizer

        @printf("\tRestart #%d: optimisation with fixed ρ=%f yields negative objective %f\n", i, makeρ(initialparamρ), intermediateresult.minimum)


        # joint optimisation of scalings α and lengthscale ρ

        local finalresult = optimize(safenegativeobj, [optparamα; initialparamρ], NelderMead(), opt)

        @printf("\tRestart #%d: final ρ=%f and negative objective %f\n", i, makeρ(finalresult.minimizer[L+1]), finalresult.minimum)

        return finalresult

    end


    allresults = [getsolution(i) for i in 1:numberofrestarts]

    result     = allresults[argmin([res.minimum for res in allresults])]

    paramopt   = result.minimizer

    @printf("\n\tOverall minimum is %f\n", result.minimum)


    #---------------------------------------------------------------------
    # instantiate learned matrix and observed variance parameter
    #---------------------------------------------------------------------

    α, ρ = unpack(paramopt)

    K = delayedCovariance(kernel, α, τ, ρ, tarray)

    KSobsB = K + Sobs + B

    makematrixsymmetric!(KSobsB)


    #---------------------------------------------------------------------
    # posterior distribution for shifts b
    #---------------------------------------------------------------------

    Σpostb = (Σb\I + Q'*((Sobs + K)\Q)) \ I

    μpostb = Σpostb * ((Q' / (Sobs + K))*Y + Σb\μb)

    posterioroffsetvector = MvNormal(μpostb, makematrixsymmetric(Σpostb))


    #---------------------------------------------------------------------
    # Functions for predicting on test data
    #---------------------------------------------------------------------


    function predictTest(ttest::Union{Array{Array{Float64, 1}, 1}, Array{T} where T<:AbstractRange{S} where S<:Real})

        Q✴  = Qmatrix(length.(ttest))

        B✴  = Q * Σb * Q✴'

        B✴✴ = Q✴ * Σb * Q✴'

        # dimensions: N × Ntest
        kB✴ = delayedCovariance(kernel, α, τ, ρ, tarray, ttest) + B✴

        # Ntest × Ntest
        cB = delayedCovariance(kernel, α, τ, ρ, ttest) + B✴✴

        # full predictive covariance
        Σpred = cB - kB✴' * (KSobsB \ kB✴)

        makematrixsymmetric!(Σpred)

        Σpred = Σpred + JITTER*I

        # predictive mean

        b̄✴ = Q✴ * μb

        μpred = kB✴' * (KSobsB \ (Y - b̄)) + b̄✴

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
    # • posterior of shift parameter
    # • scale parameters α
    # • lengthscale ρ

    result.minimum, predictTest, posterioroffsetvector, α, ρ
end
