function gpccvi(tarray, yarray, stdarray; kernel = kernel, iterations = iterations, seed = 1,  ρfixed =  ρfixed, verbose = true, S = 50)

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

    f(x) = [   softplus.(x[1:1L]); x[1L+1:2L];    transformbetween.(x[2L+1:(3L-1)],0.0,  2.0)]
    
    g(x) = [invsoftplus.(x[1:1L]); x[1L+1:2L]; invtransformbetween.(x[2L+1:(3L-1)],0.0,  2.0)]
    

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
    
    elbo = elbofy(3L-1, (3L -1)* S, helper, transform = f, invtransform = g) # take S samples per dimension/parameter

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

        local θ = rand(q)

        local α, b, τ = θ[1:L], θ[L+1:2L], [0; θ[2L+1:3L-1]]


        local K = delayedCovariance(kernel, α, τ, ρfixed, tarray)

        local KSobsB = K + Sobs

        # matrix for replicating elements

        local Q✴ = Qmatrix(length.(ttest))

        # N × Ntest
        
        local kB✴ = delayedCovariance(kernel, α, τ, ρfixed, tarray, ttest)

        # Ntest × Ntest
        
        local cB = delayedCovariance(kernel, α, τ, ρfixed, ttest)

        # full predictive covariance

        local Σpred = Symmetric(cB - kB✴' * (KSobsB \ kB✴)) + 1e-8*I

        # predictive mean

        local μpred = kB✴' * (KSobsB \ (Y - (Q*b))) + (Q✴ * b)

        # draw sample and then split it per filter

        local Ntest = length.(ttest)

        local sample = rand(MvNormal(μpred, Σpred))

        return [sample[(sum(Ntest[1:(i-1)])+1):sum(Ntest[1:i])] for i in 1:L]

    end
    
    return q, samplepredict
end
