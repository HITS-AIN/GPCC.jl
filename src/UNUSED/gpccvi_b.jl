mutable struct TrackNegativeElbo{T1,F,T2}
    elbo::T1
    bestsofar_f::F
    bestsofar_x::T2
end

function (trackelbo::TrackNegativeElbo)(x)

    local e = -trackelbo.elbo(x)

    if e < trackelbo.bestsofar_f
        trackelbo.bestsofar_f = e
        trackelbo.bestsofar_x = copy(x)
    end

    return e
end

function gpccvi(tarray, yarray, stdarray; kernel = kernel, iterations = iterations, seed = 1,  ρfixed =  ρfixed, verbose = true, S = 50, τmax = 100.0)

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

    
    
    B = 100 * diagm(map(var, yarray)) # inflated prior covariance
    
    b̄ = map(mean, yarray)             # prior mean
    
    

    #---------------------------------------------------------------------
    # Functions for constraining parameters
    #---------------------------------------------------------------------

    f(x) = [   softplus.(x[1:L]);    transformbetween.(x[1L+1:(2L-1)],0.0, τmax)]
    
    g(x) = [invsoftplus.(x[1:L]); invtransformbetween.(x[2L+1:(2L-1)],0.0, τmax)]
    

    #---------------------------------------------------------------------
    # Split parameter vector into arguments
    #---------------------------------------------------------------------
  
    function unpackθ(param)

        @assert(length(param) == 2L-1)

        local α = param[0L+1:1L]

        local τ = [0; cumsum(param[1L+1:2L-1])]

        return α, τ

    end

    function unpackψ(param)

        @assert(length(param) == 2L)

        local μb = param[0L+1:1L]

        local Cb = param[1L+1:2L]

        return μb, Cb

    end

    unpack(param) = unpackθ(param[1:2L-1]), unpackψ(param[2L:end])
 

    #---------------------------------------------------------------------
    # Define objective as marginal log-likelihood and auxiliaries
    #---------------------------------------------------------------------

    function objective(α, τ, μb, Cb)

        local Σb = Diagonal(Cb.^2)

        local K = Symmetric(delayedCovariance(kernel, α, τ, ρfixed, tarray))

        local term1 = logpdf(MvNormal(Q*μb, K + Sobs), Y) - 0.5*tr((K + Sobs)\(Q*Σb*Q'))

        local term2 = logpdf(MvNormal(b̄, B), μb) - 0.5*tr(B\Σb)

        local term3 = 0.5*sum(log.(Cb.^2)) # entropy of gaussian posterior of shift b

        term1 + term2 + term3

    end
    
    
    helper(θ, ψ) = objective(unpackθ(θ)..., unpackψ(ψ)...)
    
   
    elbo = elbofy(2L-1, 2L, (2L -1)*S, helper; transform = f, invtransform = g) # take S samples per dimension/parameter

    verbose ? display(elbo) : nothing

    testelbo = elbofy(2L-1, 2L, (2L -1)*S, helper; transform = f, invtransform = g, rg = MersenneTwister(10101)) # take S samples per dimension/parameter


    #-------------------------------------------------------
    # get point estimate to start VI
    #-------------------------------------------------------

    μ₀, ψ₀ = let

        verbose ? @printf("Initialising variational inference.\n") : nothing

        local opt = Optim.Options(show_trace = true, iterations = 100, show_every = 10)

        local aux(x) = -helper(f(x[1:2L-1]), x[2L:end])

        local out = optimize(aux, 1*randn(rg, (2L-1) + 2L), NelderMead(), opt).minimizer
     
        out[1:2L-1], out[2L:end]

    end

    Cdiag = 0.1 * ones(2L-1) 
    



    #-------------------------------------------------------
    # Report progress
    #-------------------------------------------------------

    tracknegelbo = TrackNegativeElbo(elbo, Inf, zeros(2*(2L-1)+2L))

    counter = 0

    function callback(_)
        
        counter += 1
        
        if mod(counter, 200) == 1
        
            @printf("Iter %d\t elbo is %f,\t test elbo is %f\n",counter, tracknegelbo.bestsofar_f, -testelbo(tracknegelbo.bestsofar_x))
        
        else
            
            mod(counter, 5) ==1 ? @printf("Iter %d\t elbo is %f\n",counter,tracknegelbo.bestsofar_f) : nothing

        end

        return false

    end

   

    #-------------------------------------------------------
    # optimise elbo and get optimal variational parameters
    #-------------------------------------------------------

    opt = Optim.Options(show_trace = false, extended_trace = false, iterations = iterations, show_every=1, callback = callback)
   
    paramopt = optimize(tracknegelbo, [μ₀; Cdiag; ψ₀], NelderMead(), opt).minimizer
    
    
    #----------------------------------------------
    # instantiate approximate posteriors
    #----------------------------------------------

    θ = paramopt[1:2(2L-1)]
    
    qατ = transformedgaussianposterior(elbo, θ)
    
    qb = let
        
        local ψ = paramopt[2(2L-1)+1:end]

        local μb, Cb = unpackψ(ψ)
  
         MvNormal(μb, Diagonal(Cb.^2))   

    end


    #----------------------------------------------
    # Draw samples from predictive distribution
    #----------------------------------------------
    
    samplepredict(ttest0::Array{T, 1}) where T<:Real = samplepredict([ttest0 for _ in 1:L])


    function samplepredict(ttest) 

        local θ = rand(qατ)

        local b = rand(qb)

        local α, τ = θ[1:L], [0; θ[1L+1:2L-1]]

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
    
    return qατ, qb, samplepredict
end
