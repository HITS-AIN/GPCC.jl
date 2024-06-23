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

    # get point estimate to start VI

    μ₀ = let

        verbose ? @printf("Initialising variational inference.\n") : nothing

        local opt = Optim.Options(show_trace = true, iterations = 3_000, show_every = 25)

        local aux(x) = -helper(f(x))

        optimize(aux, 1*randn(rg, 3L-1), NelderMead(), opt).minimizer
     
    end

    
    # setup progress bar 

    pr = Progress(iterations; showspeed=true, enabled = verbose)
    
    callback(st::OptimizationState) = (next!(pr; showvalues = [(:negative_elbo, st.value)]); false)
    

    # setup optimisation options for variational inference

    opt = Optim.Options(show_trace = true, iterations = iterations, show_every=1)# callback = callback)


    # initial covariance root is spherical, radius is optimised below in one-dimensional optimisation problem

    Cdiag = let

        verbose ? @printf("Initialising covariance.\n") : nothing

        local r_range = 0.1:0.1:1.0
        
        local bestindex = argmax([elbo(μ₀, r * ones(3L-1)) for r in r_range])
   
        @printf("Best r is %f\n", r_range[bestindex])

        r_range[bestindex] * ones(3L-1) 

    end

    # optimise elbo and get optimal variational parameters
   
    θ = optimize(x -> -elbo(x), [μ₀; Cdiag], NelderMead(), opt).minimizer
    

    # return approximate posterior
    
    transformedgaussianposterior(elbo, θ)
    
end
