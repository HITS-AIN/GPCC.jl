"""
    loglikel, α, postb, pred = gpcc(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, numberofrestarts = 1, initialrandom = 10, verbose = false)

Fit Gaussian Process Cross Correlation (GPCC) model for a given vector of delays `delays` and lengthscale `ρfixed``.

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
- `rng`: Set random seed generator.
- `numberofrestarts`: number of times to repeat optimisation in order to avoid suboptimal solutions due to poor initialisation (default is 1).
- `initialrandom`: number of initial solutions out of which the best is selected (default is 10).
- `ρfixed`: GP lengthscale, either set by hand or obtained by `infercommonlengthscale`.
- `verbose`: true / false (default). If set to `true`, auxiliary messages will be printed out 


Returned arguments
==================
- `loglikel`: log-likelihood reached when optimising GPCC.
- `α`: coefficients by which the latent Gaussian process is scaled in each band
- `postb`: Gaussian posterior for shift parameters returned as an object of type `MvNormal`.
- `predict`: function for predicting on out-of-sample data.

# Example
```julia-repl
julia> tobs, yobs, σobs, truedelays = simulatedata(); # produce synthetic data
julia> loglikel, α, postb, pred = gpcc(tobs, yobs, σobs; kernel = GPCC.rbf, delays = truedelays, iterations = 1000);  # fit GPCC
julia> trange = collect(-10:0.1:25); # define time interval for predictions
julia> μpred, σpred = pred(trange) # obtain predictions
julia> type(μpred), size(μpred) # predictions are also arrays of arrays, organised just like the data
julia> using PyPlot # this must be independently installed
julia> plot(trange, μpred[1], "b") # plot mean predictions for 1st band
julia> fill_between(trange, μpred[1].+σpred[1], μpred[1].-σpred[1], color="b", alpha=0.3) # plot uncertainties for 1st band
```
"""
function gpcc(tarray, yarray, stdarray; kernel = kernel, delays = delays, iterations = iterations, rng = AbstractRNG=Random.GLOBAL_RNG, numberofrestarts = 1, initialrandom = 10, verbose = false, ρfixed = ρfixed, JITTER = 1e-8)

    # Same function as below, but easier name for user to call

    gpccfixdelay(tarray, yarray, stdarray; kernel = kernel, τ = delays, iterations = iterations, rng = rng, numberofrestarts = numberofrestarts, initialrandom = initialrandom, verbose = verbose, ρfixed = ρfixed, JITTER = JITTER)


end


function gpccfixdelay(tarray, yarray, stdarray; kernel = kernel, τ = τ, iterations = iterations, rng = rng, numberofrestarts = numberofrestarts, initialrandom = initialrandom, verbose = verbose, ρfixed = ρfixed, JITTER = 1e-8)

    #---------------------------------------------------------------------
    # Check dimensions
    #---------------------------------------------------------------------

    L = length(tarray) # number of lightcurves / bands / filters

    @assert(L == length(τ) == length(yarray) == length(tarray) == length(stdarray))


    #---------------------------------------------------------------------
    # Auxiliary matrices
    #---------------------------------------------------------------------

    Y = reduce(vcat, yarray)                   # concatenated fluxes

    Q = Qmatrix(length.(tarray))               # matrix for replicating elements

    Sobs = Diagonal(reduce(vcat, stdarray).^2) # observed noise matrix


    # Prior for shift vector b

    μb = map(mean, yarray)

    Σb = 100 * Diagonal(map(var, yarray)) # inflated prior covariance

    B  = Q * Σb * Q' # prior cov after replicating to match number of observations

    b̄  = Q * μb # prior mean after replicating to match number of observations

    SobsB = Sobs + B # Additional covariance structure to be added on GP covariance


    #---------------------------------------------------------------------
    # Define objective as marginal log-likelihood and auxiliaries
    #---------------------------------------------------------------------

    delayedx = reduce(vcat, [x.-d for (x, d) in zip(tarray, τ)])

    K₁ = covariance_unit_amplitude(kernel, ρfixed, delayedx) + JITTER*I# stays fixed throughout!

   
    function objective(α) # same as one below, keep for numerical verification
        
        local A = Diagonal(Q*α)

        return logpdf(MvNormal(b̄, Symmetric(A*K₁*A + SobsB)), Y)

    end

    
    function fasterobjective(α) # same as one above, but slightly faster
        
        local A = Diagonal(Q*α)

        local C = cholesky(Symmetric(A*K₁*A + SobsB)).L

        -0.5*sum(abs2.(C\(Y-b̄))) - 0.5*2*sum(log.(diag(C))) - 0.5*log(2π)*size(C,1)

    end
      
    # Define negative objective

    negativeobjective(x) = - fasterobjective(x)


    #---------------------------------------------------------------------
    # Returns random unconstrained values for initial scaling vector α
    #---------------------------------------------------------------------

    sampleunconstrainedsolution() = map(std, yarray)  .* randn(rng, L)


    #---------------------------------------------------------------------
    # Function below calls optimiser
    #---------------------------------------------------------------------

    function getsolution()

        local opt = Optim.Options(show_trace = verbose, iterations = iterations, show_every = 2, g_tol=1e-8)

        local randomsolutions = [sampleunconstrainedsolution() for _ in 1:initialrandom]

        local bestindex = argmin(map(negativeobjective, randomsolutions))

        optimize(negativeobjective, randomsolutions[bestindex], NelderMead(), opt)

    end


    #---------------------------------------------------------------------
    # Restart optimisation multiple times as specified in `numberofrestarts`
    #---------------------------------------------------------------------

    allresults = [getsolution() for _ in 1:numberofrestarts]

    result = allresults[argmin([res.minimum for res in allresults])]

    αopt = result.minimizer

    verbose ? @printf("\n\tOverall minimum is %f\n", result.minimum) : nothing


    #---------------------------------------------------------------------
    # instantiate learned kernel matrix
    #---------------------------------------------------------------------

    A = Diagonal(Qmatrix(length.(tarray))*αopt)

    KSobsB = Symmetric(A*K₁*A + SobsB)
  

    #---------------------------------------------------------------------
    # calculate in numerical stable posterior distribution for shifts b
    #---------------------------------------------------------------------

    
    μpostb, Σpostb = let
        
        # We need to calculate  
        # Σpost = (Σb⁻¹ + Qᵀ(Sobs + AKA)⁻¹Q)⁻¹
        # and
        # μpost = Σpost{Qᵀ(Sobs + AKA)⁻¹Y + Σb⁻¹μb}

        # Naive calculation reads:
        # Symmetric((inv(Σb) + Q'*((Sobs + A*K₁*A)\Q)) \ I)

      
        # Take care of middle term, calculate its inverse using cholesky:
        local C = cholesky(Symmetric(A*K₁*A + Sobs + Q*Σb*Q')).L

        # part of woodbury
        local F = (C\Q)*Σb
        
        # this comes from woodbury
        local Σposterior = Symmetric(Σb - F'F)

        local μposterior = Σposterior*((C\Q)'*(C\Y) + Σb\μb)

        μposterior, Σposterior

    end


    #---------------------------------------------------------------------
    # Functions for predicting on test data
    #---------------------------------------------------------------------

    function predictTest(ttest::Union{Array{Array{Float64, 1}, 1}, Array{T} where T<:AbstractRange{S} where S<:Real})

        
        Q✴  = Qmatrix(length.(ttest))

        B✴  = Q * Σb * Q✴'

        B✴✴ = Q✴ * Σb * Q✴'

        # dimensions: N × Ntest
        kB✴ = delayedCovariance(kernel, αopt, τ, ρfixed, tarray, ttest) + B✴

        # Ntest × Ntest
        cB = delayedCovariance(kernel, αopt, τ, ρfixed, ttest) + B✴✴

        # full predictive covariance
        Σpred = Symmetric(cB - kB✴' * (KSobsB \ kB✴)) #+ JITTER*I

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

        Σpred = Symmetric(Σpred + Sobs✴)

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


    return -result.minimum, αopt, MvNormal(μpostb, Σpostb), predictTest
    
    
end