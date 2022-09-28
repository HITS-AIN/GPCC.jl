"""
    minopt, pred, α, postb, ρ = gpcc(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofinitialsolutions = 1, initialrandom = 5, rhomin = 0.1, rhomax = rhomax)

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
- `numberofinitialsolutions`: Number of times to repeat optimisation in order to avoid suboptimal solutions due to poor initialisation (default is 1).
- `initialrandom`: Before optimisation begins, a number of random solutions is sampled and the one with the highest likelihood becomes the starting point for the optimisation.
- `rhomin`: minimum value for lengthscale ρ of GP (default 0.1).
- `rhomax`: maximum value for lengthscale ρ of GP.


Returned arguments
==================
- `minopt`: negative log-likelihood reached when optimising GP hyperparameters.
- `predict`: function for predicting on out-of-sample data.
- `α`: coefficients by which the latent Gaussian process is scaled in each band
- `postb`: Gaussian posterior for shift parameters returned as an object of type `MvNormal`.
- `ρ`: length scale of latent Gaussian Process

# Example
```julia-repl
julia> tobs, yobs, σobs, truedelays = simulatedata(); # produce synthetic data
julia> minopt, pred, α, postb, ρ = gpcc(tobs, yobs, σobs; kernel = GPCC.matern32, delays = truedelays, iterations = 1000);  # fit GPCC
julia> trange = collect(-10:0.1:25); # define time interval for predictions
julia> μpred, σpred = pred(trange) # obtain predictions
julia> type(μpred), size(μpred) # predictions are also arrays of arrays, organised just like the data
julia> plot(trange, μpred[1], "b") # plot mean predictions for 1st band
julia> fill_between(trange, μpred[1].+σpred[1], μpred[1].-σpred[1], color="b", alpha=0.3) # plot uncertainties for 1st band
```
"""
function gpccvi(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, seed = 1, numberofinitialsolutions = 1, initialrandom = 5, rhomin = 0.1, rhomax = rhomax)

    # Same function as below, but easier name for user to call

    gpccfixdelayvi(tarray, yarray, stdarray; kernel = kernel, τ = delays, iterations = iterations, seed = seed, numberofinitialsolutions = numberofinitialsolutions, initialrandom = initialrandom, ρmin = rhomin, ρmax = rhomax)


end


function gpccfixdelayvi(tarray, yarray, stdarray; kernel = kernel, τ = τ, iterations = iterations, seed = 1, numberofinitialsolutions = 1, initialrandom = 5, ρmin = 0.1, ρmax = 20.0)

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


    μb = map(mean, yarray)             # prior mean

    Σb = 100 * diagm(map(var, yarray)) # inflated prior covariance

    B  = Q * Σb * Q'

    b̄  = Q * μb

    #---------------------------------------------------------------------
    # Let user know what is being run
    #---------------------------------------------------------------------

    # informuser(seed = seed, iterations = iterations, numberofinitialsolutions = numberofinitialsolutions,
    #             initialrandom = initialrandom, JITTER = JITTER, ρmin = ρmin, ρmax = ρmax, Σb = Σb)


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

        local K = delayedCovariance(kernel, α, τ, ρ, tarray) + Sobs + B

        makematrixsymmetric!(K)

        return logpdf(MvNormal(b̄, K), Y)

    end

    # convenient call

    objective(param) = objective(unpack(param)...)

    # Define negative objective

    negativeobjective(x) = - objective(x)

    # Auxiliary objective catches exceptions

    safenegativeobj = safewrapper(negativeobjective)

    # safeobj = safewrapper(objective)


    #---------------------------------------------------------------------
    # Define initial values for lengthscale ρ
    #---------------------------------------------------------------------

    initialρvalues = let

        if numberofinitialsolutions == 1 || numberofinitialsolutions == 2

            # pick initial ρ values randomly

            rand(rg, Uniform(ρmin + 1e-3, ρmax - 1e-3), numberofinitialsolutions)

        else

            # initial ρ values on grid

            collect(MiscUtil.logrange(ρmin + 1e-3, ρmax - 1e-3, numberofinitialsolutions))

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
    # Evaluate initial solutions
    #---------------------------------------------------------------------

    function getinitialsolution()

        local randomsolutions = [sampleunconstrainedsolution(i) for i in 1:numberofinitialsolutions]

        local bestindex = argmin(map(safenegativeobj, randomsolutions))

        return randomsolutions[bestindex]
        
    end


    #---------------------------------------------------------------------
    # 
    #---------------------------------------------------------------------

    
    paramopt0   = getinitialsolution()

    
    #---------------------------------------------------------------------
    # 
    #---------------------------------------------------------------------

    post, mll = ApproximateVI.VI(objective, paramopt0, optimiser = Optim.NelderMead(), iterations = iterations, show_every=50, S=75)

    

    #---------------------------------------------------------------------
    # instantiate learned kernel matrix
    #---------------------------------------------------------------------




    #---------------------------------------------------------------------
    # posterior distribution for shifts b
    #---------------------------------------------------------------------

    # Σpostb = (Σb\I + Q'*((Sobs + K)\Q)) \ I

    # μpostb = Σpostb * ((Q' / (Sobs + K))*Y + Σb\μb)

    # postb = MvNormal(μpostb, makematrixsymmetric(Σpostb))


    #---------------------------------------------------------------------
    # Functions for predicting on test data
    #---------------------------------------------------------------------

    function predictTest(ttest::Union{Array{Array{Float64, 1}, 1}, Array{T} where T<:AbstractRange{S} where S<:Real})


        paramopt = rand(post)

        α, ρ = unpack(paramopt)

        K = delayedCovariance(kernel, α, τ, ρ, tarray)
    
        KSobsB = K + Sobs + B
    
        makematrixsymmetric!(KSobsB)



        
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
    # • optimised free parameters

    mll, predictTest, post, unpack
end
