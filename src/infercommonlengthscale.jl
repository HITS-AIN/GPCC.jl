function infercommonlengthscale(tarray, yarray, stdarray; kernel = kernel, iterations = iterations, seed = 1, numberofrestarts = 1, initialrandom = 5, ρmin = 0.1, ρmax = 20.0, verbose = verbose)

    #---------------------------------------------------------------------
    # Fix random seed for reproducibility
    #---------------------------------------------------------------------

    rg = MersenneTwister(seed)


    #---------------------------------------------------------------------
    # Set constants
    #---------------------------------------------------------------------

    JITTER = 1e-8

    L = length(tarray)


    Sobs = [Diagonal(stdarray[l].^2) for l in 1:L] # observed noise matrix


    #---------------------------------------------------------------------
    # Check dimensions
    #---------------------------------------------------------------------

    @assert(L == length(yarray) == length(tarray) == length(stdarray))


    covmatrix(x, y, α, ρ) = [α*α * kernel(x₁, x₂ ; ρ=ρ)  for x₁ in x, x₂ in y]


    #---------------------------------------------------------------------
    # Functions for constraining parameters
    #---------------------------------------------------------------------

    makeα(x) = softplus(x)

    makeρ(x) = transformbetween(x, ρmin, ρmax)

    function unpack(param)

        @assert(length(param) == 2L + 1)

        local α = makeα.(param[0L+1:1L])

        local b = param[1L+1:2L]

        local ρ = makeρ(param[2L+1])

        return α, b, ρ

    end


    #---------------------------------------------------------------------
    # Define objective as marginal log-likelihood and auxiliaries
    #---------------------------------------------------------------------

    function objective(α, b, ρ)

        local aux = zero(eltype(α))

        for l in 1:L

            local K = Symmetric(covmatrix(tarray[l], tarray[l], α[l], ρ) + Sobs[l])

            aux += logpdf(MvNormal(ones(length(yarray[l])) * b[l], K), yarray[l])

        end

        return aux

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


    if verbose 
        
        @printf("\n\tInitial ρ values are:\n")

        map(x -> @printf("\t%f\n", x), initialρvalues)

    end


    #---------------------------------------------------------------------
    # Returns random values for initial scaling vector α and shift vector v
    #---------------------------------------------------------------------

    sampleα() = map(var, yarray)  .* (rand(rg, L) * (1.2 - 0.8) .+ 0.8)

    sampleb() = map(mean, yarray) .* (rand(rg, L) * (1.2 - 0.8) .+ 0.8)


    #---------------------------------------------------------------------
    # Returns random unconstrained solution
    #---------------------------------------------------------------------

    sampleunconstrainedsolution(i) = [invmakepositive.(sampleα()); sampleb();
                                      invtransformbetween(initialρvalues[i], ρmin, ρmax)]


    #---------------------------------------------------------------------
    # Function below calls optimiser
    #---------------------------------------------------------------------

    function getsolution(i)

        local opt = Optim.Options(show_trace = verbose, iterations = iterations, show_every = 2, g_tol=1e-6)

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

    if verbose
        @printf("\n\tOverall minimum is %f\n", result.minimum)
    end


    #---------------------------------------------------------------------
    # instantiate learned kernel matrix
    #---------------------------------------------------------------------

    α, b, ρ = unpack(paramopt)

    K = [covmatrix(tarray[l], tarray[l], α[l], ρ) + Sobs[l] for l in 1:L]

    #---------------------------------------------------------------------
    # Functions for predicting on test data
    #---------------------------------------------------------------------

    function predictTest(l, ttest)

        # dimensions: N × Ntest
        # kB✴ = delayedCovariance(kernel, α, τ, ρ, tarray, ttest)
        kB✴ = covmatrix(tarray[l], ttest, α[l], ρ)

        # Ntest × Ntest
        # cB = delayedCovariance(kernel, α, τ, ρ, ttest) + B✴✴
        cB = covmatrix(ttest, ttest, α[l], ρ)

        # full predictive covariance
        Σpred = Symmetric(cB - kB✴' * (K[l] \ kB✴) + JITTER*I)

        
        μpred = kB✴' * (K[l] \ (yarray[l] .- b[l])) .+ b[l]

        return μpred, Σpred

    end




    # return:
    # • function value returned from optimisation
    # • prediction function
    # • optimised free parameters

    -result.minimum, predictTest, (α, b, ρ)
end
