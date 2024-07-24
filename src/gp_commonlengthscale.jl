function infercommonlengthscale(tarray, yarray, stdarray; kernel = kernel, iterations = iterations, rng = AbstractRNG=Random.GLOBAL_RNG, numberofrestarts = 1, initialrandom = 5, ρmin = 0.1, ρmax = 20.0, verbose = verbose, JITTER = 1e-8)

    gp_commonlengthscale(tarray, yarray, stdarray; kernel = kernel, iterations = iterations, rng = rng, numberofrestarts = numberofrestarts, initialrandom = initialrandom, ρmin = ρmin, ρmax = ρmax, verbose = verbose, JITTER = JITTER)[4]


end

function gp_commonlengthscale(tarray, yarray, stdarray; kernel = kernel, iterations = iterations, rng = AbstractRNG=Random.GLOBAL_RNG, numberofrestarts = 1, initialrandom = 5, ρmin = 0.1, ρmax = 20.0, verbose = false, JITTER = 1e-8)

    #---------------------------------------------------------------------
    # Set constants
    #---------------------------------------------------------------------

    L = length(tarray)

    Sobs = [Diagonal(stdarray[l].^2) for l in 1:L] # observed noise matrix


    # Prior for shift vector b

    𝟏 = [ones(length(tarray[l])) for l in 1:L]

    μb  = map(mean, yarray)

    σ²b = 100 * Diagonal(map(var, yarray)) # inflated prior variance

    B  = [σ²b[l]*𝟏[l]*𝟏[l]' for l in 1:L] # prior cov after replicating to match number of observations

    # SobsB = [Sobs[l] + B[l] for l in 1:L] # combined covariance to be added to GP covariance


    #---------------------------------------------------------------------
    # Check dimensions
    #---------------------------------------------------------------------

    @assert(L == length(yarray) == length(tarray) == length(stdarray))


    covmatrix(x, y, α, ρ) = [α * α * kernel(x₁, x₂ ; ρ=ρ)  for x₁ in x, x₂ in y]


    #---------------------------------------------------------------------
    # Functions for constraining parameters
    #---------------------------------------------------------------------

    makeρ(x) = transformbetween(x, ρmin, ρmax)

    function unpack(param)

        @assert(length(param) == 1L + 1)

        local α = param[0L+1:1L]

        local ρ = makeρ(param[1L+1])

        return α, ρ

    end


    #---------------------------------------------------------------------
    # Define objective as marginal log-likelihood and auxiliaries
    #---------------------------------------------------------------------

    function objective(α, ρ)

        local aux = zero(eltype(α))

        for l in 1:L

            local K = Symmetric(covmatrix(tarray[l], tarray[l], α[l], ρ) + Sobs[l] + B[l])

            aux += logpdf(MvNormal(μb[l]*𝟏[l], K), yarray[l])

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

            rand(rng, Uniform(ρmin + 1e-3, ρmax - 1e-3), numberofrestarts)

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

    sampleα() = map(std, yarray)  .* randn(rng, L)


    #---------------------------------------------------------------------
    # Returns random unconstrained solution
    #---------------------------------------------------------------------

    sampleunconstrainedsolution(i) = [sampleα();
                                      invtransformbetween(initialρvalues[i], ρmin, ρmax)]


    #---------------------------------------------------------------------
    # Function below calls optimiser
    #---------------------------------------------------------------------

    function getsolution(i)

        local opt = Optim.Options(show_trace = verbose, iterations = iterations, show_every = 100, g_tol=1e-6)

        local randomsolutions = [sampleunconstrainedsolution(i) for _ in 1:initialrandom]

        local bestindex = argmin(map(safenegativeobj, randomsolutions))

        local finalresult = optimize(safenegativeobj,randomsolutions[bestindex], NelderMead(), opt)

        return finalresult

    end


    #---------------------------------------------------------------------
    # Restart optimisation multiple times as specified in `numberofrestarts`
    #---------------------------------------------------------------------

    allresults = @showprogress "optimising length scale" [getsolution(i) for i in 1:numberofrestarts]

    result     = allresults[argmin([res.minimum for res in allresults])]

    paramopt   = result.minimizer

    if verbose
        @printf("\n\tOverall minimum is %f\n", result.minimum)
    end


    #---------------------------------------------------------------------
    # instantiate learned kernel matrix
    #---------------------------------------------------------------------

    α, ρ = unpack(paramopt)

    K = [covmatrix(tarray[l], tarray[l], α[l], ρ) + Sobs[l] + B[l] for l in 1:L]

    #---------------------------------------------------------------------
    # Functions for predicting on test data
    #---------------------------------------------------------------------

    function predictTest(l, ttest)

        𝟏✴ = ones(length(ttest))

        B✴  = σ²b[l] * 𝟏[l]*𝟏✴'

        B✴✴ = σ²b[l] * 𝟏✴*𝟏✴'

        # dimensions: N × Ntest
        K✴ = covmatrix(tarray[l], ttest, α[l], ρ) + B✴

        # Ntest × Ntest
        K✴✴  = covmatrix(ttest, ttest, α[l], ρ) + B✴✴

        # full predictive covariance - see 2.26 in RW
        Σpred = Symmetric(K✴✴ - K✴' * (K[l] \ K✴))

        # predictive mean - see 2.25, 2.38 and 2.41 in RW
        μpred = K✴' * (K[l] \ (yarray[l] .- μb[l])) .+ μb[l]

        return μpred, Σpred

    end




    # return:
    # • function value returned from optimisation
    # • prediction function
    # • optimised scaling parameters
    # • optimised length scale

    -result.minimum, predictTest, α, ρ

end
