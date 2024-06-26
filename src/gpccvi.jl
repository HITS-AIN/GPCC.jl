function gpccvi(tarray, yarray, stdarray; kernel = kernel, iterations = iterations, seed = 1,  ρfixed =  ρfixed, verbose = true)

    #---------------------------------------------------------------------
    # Fix random seed, get number of filters and check dimesions
    #---------------------------------------------------------------------

    rg = MersenneTwister(seed)
    
    L = length(tarray)
    
    @assert(L == length(yarray) == length(tarray) == length(stdarray))


    #---------------------------------------------------------------------
    # Auxiliary matrices
    #---------------------------------------------------------------------

    Y = reduce(vcat, yarray)                   # concatenated fluxes

    Q = Qmatrix(length.(tarray))               # matrix for replicating elements

    Sobs = Diagonal(reduce(vcat, stdarray).^2) # observed noise matrix

    
    #---------------------------------------------------------------------
    # Functions for constraining parameters
    #---------------------------------------------------------------------

    f(x) = [   softplus.(x[1:1L]); x[1L+1:2L];    softplus.(x[2L+1:(3L-1)])]
    
    g(x) = [invsoftplus.(x[1:1L]); x[1L+1:2L]; invsoftplus.(x[2L+1:(3L-1)])]
    

    #---------------------------------------------------------------------
    # Split parameter vector into arguments
    #---------------------------------------------------------------------
  
    function unpack(param)

        @assert(length(param) == 3L-1)

        local α = param[0L+1:1L]

        local b = param[1L+1:2L]

        local τ = [0; cumsum(param[2L+1:3L-1])]

        # local ρ = param[3L]

        return α, b, τ, ρfixed

    end

 
    #---------------------------------------------------------------------
    # Define objective as marginal log-likelihood and auxiliaries
    #---------------------------------------------------------------------

    function objective(α, b, τ, ρ)

        local K = Symmetric(delayedCovariance(kernel, α, τ, ρ, tarray) + Sobs)

        return logpdf(MvNormal(Q*b, K), Y)

    end
    
    helper(p) = objective(unpack(p)...)
    
    elbo = elbofy(3L-1, (3L -1)* 50, helper, transform = f, invtransform = g) # take 50 samples per dimension/parameter

    verbose ? display(elbo) : nothing


    #-------------------------------------------------------
    # get point estimate to start VI
    #-------------------------------------------------------

    μ₀ = let

        verbose ? @printf("Initialising variational inference.\n") : nothing

        local opt = Optim.Options(show_trace = true, iterations = 3_000, show_every = 25)

        local aux(x) = -helper(f(x))

        optimize(aux, 1*randn(rg, 3L-1), NelderMead(), opt).minimizer
     
    end

    
    # # setup progress bar 

    # pr = Progress(iterations; showspeed=true, enabled = verbose)
    
    # callback(st::OptimizationState) = (next!(pr; showvalues = [(:negative_elbo, st.value)]); false)
    

    
    #-------------------------------------------------------
    # initial covariance root is spherical,
    # radius is optimised below in one-dimensional optimisation problem
    #-------------------------------------------------------
    
    Cdiag = let
        
        verbose ? @printf("Initialising covariance.\n") : nothing
        
        local r_range = 0.1:0.1:1.0
        
        local bestindex = argmax([elbo(μ₀, r * ones(3L-1)) for r in r_range])
        
        @printf("Best r is %f\n", r_range[bestindex])
        
        r_range[bestindex] * ones(3L-1) 
        
    end
    
    
    #-------------------------------------------------------
    # optimise elbo and get optimal variational parameters
    #-------------------------------------------------------

    opt = Optim.Options(show_trace = true, iterations = iterations, show_every=1)# callback = callback)
   
    θ = optimize(x -> -elbo(x), [μ₀; Cdiag], NelderMead(), opt).minimizer
    

    #----------------------------------------------
    # instantiate approximate posterior
    #----------------------------------------------

    q = transformedgaussianposterior(elbo, θ)


    #----------------------------------------------
    # Draw samples from predictive distribution
    #----------------------------------------------
    
    samplepredict(ttest0::Array{T, 1}) where T<:Real = samplepredict([ttest0 for _ in 1:L])


    function samplepredict(ttest) 

        local JITTER = 1e-10

        local ρ = ρfixed

        local θ = rand(q)

        local α, b, τ = θ[1:L], θ[L+1:2L], [0; θ[2L+1:3L-1]]


        local K = delayedCovariance(kernel, α, τ, ρ, tarray)

        local KSobsB = K + Sobs

        # matrix for replicating elements
        local Q✴ = Qmatrix(length.(ttest))

        # dimensions: N × Ntest
        local kB✴ = delayedCovariance(kernel, α, τ, ρ, tarray, ttest)

        # Ntest × Ntest
        local cB = delayedCovariance(kernel, α, τ, ρ, ttest)

        # full predictive covariance
        local Σpred = Symmetric(cB - kB✴' * (KSobsB \ kB✴)) + JITTER*I

        # predictive mean

        local μpred = kB✴' * (KSobsB \ (Y - (Q*b))) + (Q✴ * b)

        Ntest = length.(ttest)
    
        # split predictions for each filter
        local μ = [μpred[(sum(Ntest[1:(i-1)])+1):sum(Ntest[1:i])] for i in 1:L]
        local Σ = [Σpred[(sum(Ntest[1:(i-1)])+1):sum(Ntest[1:i])] for i in 1:L]
        
        return [rand(MvNormal(μᵢ, Σᵢ)) for (μᵢ, Σᵢ) in zip(μ, Σ)]

    end
    
    return q, samplepredict
end
