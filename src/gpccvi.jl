"""
    gpcc(tarray, yarray, stdarray;  minimumSqLen = 0.1, iterations = 100, numberofrestarts = 1, seed = 1)

Gaussian Process Cross Correlation

The first 3 arguments to the function concern the data.
The rest of the arguments control the algorithm.


# Data arguments

`tarray` is an array of arrays.
The i-th inner array holds the observation times for the i-th band.

`yarray` is an array of arrays.
The i-th inner array holds the observed fluxes for the i-th band.

`stdarray` is an array of arrays.
The i-th inner array holds the observed standard deviations for the i-th band.


# Algorithm arguments

`minimumSqLen` is the minimum squared length scale  (ℓ²) value for RBF kernel function.

`iterations` number of iterations for optimiser to run.

`numberofrestarts` number of times that the optimisation is re-run in order to avoid local minima.

`seed` controls the random number generator which us used to sample the initial solutions randomly.

# Output



# Examples
```julia-repl
julia> using PyPlot
julia> t, y, σ = simulatedata() # generate synthetic data
julia> pred = gpcc(t, y, σ; minimumSqLen=0.1, numberofrestarts = 10, iterations=10_000)
julia> xtest = -13:0.1:26 # define an interval on which we predict
julia> μ₁, μ₂ = pred(xtest)[1] # get mean predictions
julia> σ₁, σ₂ = pred(xtest)[2] # get std predictions
julia> plot(xtest, μ₁, "k-",  label="mean pred for 1st time series")
julia> plot(xtest, μ₂, "k--", label="mean pred for 2nd time series")
julia> fill_between(xtest, μ₁.+σ₁, μ₁.-σ₁, color="k", alpha=0.1)
julia> fill_between(xtest, μ₂.+σ₂, μ₂.-σ₂, color="k", alpha=0.1)
julia> legend()
```
"""
function gpccvi(tarray, yarray, stdarray; minimumSqLen = 0.1, numberofrestarts = 1, iterations = 100, seed = 1)

    rg = MersenneTwister(seed)

    JITTER = 1e-6

    minimumScale = 0.1


    #---------------------------------------------------------------------
    # Sort out dimensions
    #---------------------------------------------------------------------

    N = map(length, tarray)

    L = length(yarray); @assert(L == length(yarray) == length(tarray) == length(stdarray))

    Y = reduce(vcat, yarray)

    Σobs = Diagonal(reduce(vcat, stdarray).^2)


    #---------------------------------------------------------------------
    # Report
    #---------------------------------------------------------------------

    # let user know what is run
    str = "\nRunning gpccvi\n\n"
    print(Crayon(foreground = :light_magenta, bold = true), @sprintf("%s", str), Crayon(reset = true))

    @printf("\t \t number of restarts is %d\n", numberofrestarts)
    @printf("\t \t iterations are %d\n", iterations)
    @printf("\t \t number of bands is %d\n", L)
    @printf("\t \t jitter is %e\n", JITTER)
    @printf("\t \t random seed is %e (UNUSED!)\n", seed)
    @printf("\t \t minimum scale is %f\n", minimumScale)
    @printf("\t \t minimum ℓ² is %f\n\n", minimumSqLen)


    #---------------------------------------------------------------------
    function unpack(param)
    #---------------------------------------------------------------------

        @assert(length(param) == 3L)

        local mark = 0

        local delays = [0.0; (param[mark + 1:mark + L-1])]

        mark += L-1

        local scale  = exp.(param[mark + 1:mark + L]) .+ minimumScale

        mark += L

        local b      = param[mark + 1:mark + L]

        mark += L

        local ℓ²     = exp(param[mark + 1]) + minimumSqLen

        mark += 1

        @assert(mark == 3L)

        return delays, scale, b, ℓ²

    end


    #---------------------------------------------------------------------
    function objective(param)
    #---------------------------------------------------------------------

        local delays, scale, b, ℓ² = unpack(param)

        local C = delayedCovariance(scale, delays, ℓ², tarray) + Σobs + JITTER*I

        makematrixsymmetric!(C)

        local mu = reduce(vcat, [ones(N[i])*b[i] for i in 1:L])

        return logpdf(MvNormal(mu, C), Y)

    end


    #---------------------------------------------------------------------
    # Objective wrapped in function that catches posdef exception
    #---------------------------------------------------------------------

    safeobj = safewrapper(objective)


    #---------------------------------------------------------------------
    # Call optimiser
    #---------------------------------------------------------------------


    getsolution(x) = optimize(x->-safeobj(x), x, NelderMead(), Optim.Options(show_trace = true, iterations = 10_000, show_every = 100, g_tol=1e-6), autodiff=:forward)

    samplescale() = rand(rg)*(5-minimumScale) + minimumScale

    sampleℓ²()    = rand(rg)*(30-minimumSqLen) + minimumSqLen

    initialrandomsolutions = [[log.(rand(rg, L-1)*5 .+ 0.1);
                                [log(samplescale()) for l in 1:L];
                                randn(rg, L)*10;
                                log(sampleℓ²()-minimumSqLen)] for i in 1:numberofrestarts]


    results = map(getsolution, initialrandomsolutions)

    bestindex = argmin([r.minimum for r in results])

    result = results[bestindex]

    bestinitialsolution = result.minimizer

    post, mll = MVI(safeobj, bestinitialsolution, optimiser = Optim.LBFGS(), iterations = iterations, show_every=1, S=200)

    return post, mll, unpack
    #---------------------------------------------------------------------
    # instantiate learned matrix and observed variance parameter
    #---------------------------------------------------------------------

    @show delays, scale, b, ℓ² = unpack(paramopt)

    mu = reduce(vcat, [ones(N[i])*b[i] for i in 1:L])

    C = delayedCovariance(scale, delays, ℓ², tarray)

    Kσ2 = C + Σobs + JITTER*I

    makematrixsymmetric!(Kσ2)


    #---------------------------------------------------------------------
    # Functions for predicting on test data
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------------
    function predictTest(ttest::Array{Array{Float64, 1}, 1})
    #---------------------------------------------------------------------------

        local Ntest = length(ttest)

        # dimensions: N × Ntest
        # k = bigcovariancemodified(garray, scale, ℓ², tarray, ttest)
        local k = delayedCovariance(scale, delays, ℓ², tarray, ttest)

        # Ntest × 1
        # c = bigcovariancemodified(garray, scale, ℓ², ttest)
        local c = delayedCovariance(scale, delays, ℓ², ttest)

        # full predictive covariance and mean
        local Σpred = c-k'*(Kσ2\k) + JITTER*I

        makematrixsymmetric!(Σpred)


        local mupred = reduce(vcat, [ones(length(ttest[i]))*b[i] for i in 1:L])


        local μpred = k'*(Kσ2\(Y.-mu))

        @assert(length(mupred) == length(μpred))

        μpred = μpred .+ mupred


        return μpred, Σpred

    end


    #---------------------------------------------------------------------------
    function predictTest(ttest::Union{AbstractRange{Float64}, Array{Float64,1}})
    #---------------------------------------------------------------------------

        local Ntest = length(ttest)

        # dimensions: N × Ntest

        # k = bigcovariancemodified(garray, scale, ℓ², tarray, [ttest for _ in 1:L])
        local k = delayedCovariance(scale, delays, ℓ², tarray, [ttest for _ in 1:L])

        # Ntest × 1
        # c = bigcovariancemodified(garray, scale, ℓ², [ttest for _ in 1:L])
        local c = delayedCovariance(scale, delays, ℓ², [ttest for _ in 1:L])

        # full predictive covariance and mean

        local Σpred = c-k'*(Kσ2\k) + JITTER*I

        makematrixsymmetric!(Σpred)


        local μpred = k'*(Kσ2\(Y.-mu))

        # return predictions per "band" and collapse full covariance to diagonal only

        local μ_per_band = [μpred[idx] for idx in Iterators.partition(1:L*Ntest, Ntest)]

        local σ_per_band = [sqrt.(max.(diag(Σpred)[idx], 1e-6)) for idx in Iterators.partition(1:L*Ntest, Ntest)]

        for i in 1:length(μ_per_band)

            μ_per_band[i] .+= b[i]

        end

        return μ_per_band, σ_per_band

    end


    #---------------------------------------------------------------------------
    function predictTest(ttest::Array{Array{Float64, 1}, 1},
                         ytest::Array{Array{Float64, 1}, 1},
                         σtest::Array{Array{Float64, 1}, 1})
    #---------------------------------------------------------------------------

        local μpred, Σpred = predictTest(ttest)

        local Σobs = Diagonal(reduce(vcat, σtest).^2)

        Σpred = Σpred + Σobs + JITTER*I

        makematrixsymmetric!(Σpred)

        try

            return logpdf(MvNormal(μpred, Σpred), reduce(vcat, ytest))

        catch exception

            if isa(exception, PosDefException)

                local Evalues, Evector = eigen(Σpred)

                local minimumeigenvalue = 1e-6

                local newΣpred = Evector * Diagonal(max.(Evalues, minimumeigenvalue)) * Evector'

                makematrixsymmetric!(newΣpred)

                return logpdf(MvNormal(μpred, newΣpred), reduce(vcat, ytest))

            else

                throw(exception)

            end

        end

    end


    # return prediction functions
    predictTest

end
