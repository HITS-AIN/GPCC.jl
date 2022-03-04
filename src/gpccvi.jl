function gpcc2vi(tarray, yarray, stdarray; iterations = iterations, seed = 1, numberofrestarts = 1, ρmin = 0.1)

    #---------------------------------------------------------------------
    # Fix random seed for reproducibility
    #---------------------------------------------------------------------

    rg = MersenneTwister(seed)


    #---------------------------------------------------------------------
    # Set constants
    #---------------------------------------------------------------------

    JITTER = 1e-8   # JITTER to avoid numerical instability

    Tmax = maximum(reduce(vcat, tarray))

    L = length(tarray)


    #---------------------------------------------------------------------
    # Sort out dimensions
    #---------------------------------------------------------------------

    Narray = length.(tarray)

    @assert(L == length(yarray) == length(tarray) == length(stdarray))


    #---------------------------------------------------------------------
    # Auxiliary matrices
    #---------------------------------------------------------------------

    Y = reduce(vcat, yarray)          # concatenated time series observations

    Sobs = Diagonal(reduce(vcat, stdarray).^2)

    Q  = Qmatrix(Narray)              # matrix for replicating elements

    μb = map(mean, yarray)            # prior mean

    Σb = 100 * diagm(map(var, yarray)) # inflated prior covariance

    B  = Q * Σb * Q'

    b̄  = Q * μb


    #---------------------------------------------------------------------
    # Let user know what is being run
    #---------------------------------------------------------------------

    informuser(seed = seed, iterations = iterations, numberofrestarts = numberofrestarts,
                JITTER = JITTER, ρmin = ρmin, Σb = Σb)


    #---------------------------------------------------------------------
    function unpack(param)
    #---------------------------------------------------------------------

        @assert(length(param) == 2L)

        local scale  = makepositive.(param[1+0L:1*L]) .+ 1e-3

        local delays = [0.0; makepositive.(param[1+1L:2*L-1])]

        local ρ     = makepositive(param[2L]) + 1e-3

        return scale, delays, ρ

    end


    #---------------------------------------------------------------------
    function objective(param)
    #---------------------------------------------------------------------

        local scale, delays, ρ = unpack(param)

        local K = delayedCovariance(scale, delays, ρ, tarray)

        local KSobsB = K + Sobs + B

        makematrixsymmetric!(KSobsB)

        return logpdf(MvNormal(b̄, KSobsB), Y)

    end

    negativeobjective(x) = - objective(x)

    safenegativeobj = safewrapper(negativeobjective)

    safeobj = safewrapper(objective)

    #---------------------------------------------------------------------
    # Call optimiser and initialise with random search
    #---------------------------------------------------------------------

    function getsolution()

        sampleρ()       = rand(rg, Uniform(ρmin, Tmax))

        samplescales()  = map(var, yarray) .* (rand(rg, L) * (1.2 - 0.8) .+ 0.8)

        sampledelays()  = cumsum(rand(rg, L-1) * Tmax)

        randomsolutions = [[invmakepositive.(samplescales());
                            invmakepositive.(sampledelays());
                            invmakepositive(sampleρ())] for i in 1:50]

        bestindex = argmin(map(safenegativeobj, randomsolutions))

        opt = Optim.Options(show_trace = true, iterations = iterations, show_every = 10, g_tol=1e-8)

        optimize(safenegativeobj, randomsolutions[bestindex], NelderMead(), opt)

    end

    allresults = [getsolution() for _ in 1:numberofrestarts]

    result = allresults[argmin([res.minimum for res in allresults])]

    paramopt = result.minimizer

    post, mll = VI(safeobj, paramopt, optimiser = Optim.LBFGS(), iterations = iterations, show_every=1, S=75)

    return mll, post, unpack

    #---------------------------------------------------------------------
    # instantiate learned matrix and observed variance parameter
    #---------------------------------------------------------------------

    @show scale, delays, ρ = unpack(paramopt)

    K = delayedCovariance(scale, delays, ρ, tarray)

    KSobsB = K + Sobs + B

    makematrixsymmetric!(KSobsB)


    #---------------------------------------------------------------------
    # Functions for predicting on test data
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------------
    function predictTest(ttest::Union{Array{Array{Float64, 1}, 1}, Array{T} where T<:AbstractRange{S} where S<:Real})
    #---------------------------------------------------------------------------

        Q✴  = Qmatrix(length.(ttest))

        B✴  = Q * Σb * Q✴'

        B✴✴ = Q✴ * Σb * Q✴'

        # dimensions: N × Ntest
        kB✴ = delayedCovariance(scale, delays, ρ, tarray, ttest) + B✴

        # Ntest × Ntest
        cB = delayedCovariance(scale, delays, ρ, ttest) + B✴✴

        # full predictive covariance
        Σpred = cB - kB✴' * (KSobsB \ kB✴)

        makematrixsymmetric!(Σpred)

        Σpred = Σpred + JITTER*I

        # predictive mean

        b̄✴ = Q✴ * μb

        μpred = kB✴' * (KSobsB \ (Y - b̄)) + b̄✴

        return μpred, Σpred

    end


    #---------------------------------------------------------------------------
    function predictTest(ttest::Union{AbstractRange{Float64}, Array{Float64,1}})
    #---------------------------------------------------------------------------

        Ntest = length(ttest)

        local μpred, Σpred = predictTest([ttest for _ in 1:L])

        # return predictions per "band" and collapse full covariance to diagonal of standard deviations only

        μ_per_band = [μpred[idx] for idx in Iterators.partition(1:L*Ntest, Ntest)]

        σ_per_band = [sqrt.(max.(diag(Σpred)[idx], 1e-6)) for idx in Iterators.partition(1:L*Ntest, Ntest)]

        return μ_per_band, σ_per_band

    end


    #---------------------------------------------------------------------------
    function predictTest(ttest::Array{Array{Float64, 1}, 1},
                         ytest::Array{Array{Float64, 1}, 1},
                         σtest::Array{Array{Float64, 1}, 1})
    #---------------------------------------------------------------------------

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


    # return function value returned from optimisation and prediction function

    result.minimum, predictTest

end





function informuser(; seed = seed, iterations = iterations, numberofrestarts = numberofrestarts,
                    JITTER = JITTER, ρmin = ρmin, Σb = Σb)

    colourprint(@sprintf("Running gpcc2vi with random seed %d\n", seed), foreground = :light_blue, bold = true)
    @printf("\t iterations             = %d\n", iterations)
    @printf("\t numberofrestarts       = %d\n", numberofrestarts)
    @printf("\t JITTER                 = %e\n", JITTER)
    @printf("\t ρmin                   = %f\n", ρmin)
    @printf("\t ρmax                   = %f\n", ρmax)
    @printf("\t Σb                     = "); map(x->@printf("%.3f ",x), diag(Σb)); @printf("\n")
end
