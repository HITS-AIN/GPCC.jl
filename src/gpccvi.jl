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


function gpccvi(tarray, yarray, stdarray; delays = delays, kernel = kernel, iterations = iterations, seed = 1,  ρfixed =  ρfixed, verbose = true, S = 50, τmax = 100.0)

    #---------------------------------------------------------------------
    # Fix random seed, get number of filters and check dimesions
    #---------------------------------------------------------------------

    rg = MersenneTwister(seed)
    
    L = length(tarray)
    
    @assert(L == length(yarray) == length(tarray) == length(stdarray) == length(delays))


    #---------------------------------------------------------------------
    # Auxiliary matrices
    #---------------------------------------------------------------------

    Y = reduce(vcat, yarray)                   # concatenated fluxes

    Q = Qmatrix(length.(tarray))               # matrix for replicating elements

    Sobs = Diagonal(reduce(vcat, stdarray).^2) # observed noise matrix

    μb = map(mean, yarray)             # prior mean

    Σb = 10 * Diagonal(map(var, yarray)) # inflated prior covariance

    B  = Q * Σb * Q'

    b̄  = Q * μb

    SobsB = Sobs + B

    #---------------------------------------------------------------------
    # Functions for constraining parameters
    #---------------------------------------------------------------------

    f(x) = softplus.(x[1:1L])
    
    g(x) = invsoftplus.(x[1:1L])
    

    #---------------------------------------------------------------------
    # Split parameter vector into arguments
    #---------------------------------------------------------------------
  
    function unpack(param)

        @assert(length(param) == 1L)

        local α = param[0L+1:1L]

        return α

    end

 
    #---------------------------------------------------------------------
    # Define objective as marginal log-likelihood and auxiliaries
    #---------------------------------------------------------------------


    delayedx = reduce(vcat, [x.-d for (x, d) in zip(tarray, delays)])

    K₁ = covariance_unit_amplitude(kernel, ρfixed, delayedx) # stays fixed throughout!

    Y_minus_b̄ = Y-b̄

    function objective(α)
        
        local A = Diagonal(Q*α)
        
        local C = cholesky(Symmetric(A*K₁*A + SobsB)).L

        -0.5*sum(abs2.(C\(Y_minus_b̄))) - 0.5*2*sum(log.(diag(C))) #- 0.5*log(2π)*size(C,1)

    end
    
    
    helper(p) = objective(unpack(p))
    
    elbo = elbofy(1*L, 1*L*S, helper, transform = f, invtransform = g) # take S samples per dimension/parameter


    testelbo = elbofy(1*L, 1*L*S, helper; transform = f, invtransform = g, rg = MersenneTwister(10101)) # take S samples per dimension/parameter

    verbose ? display(elbo) : nothing


    #-------------------------------------------------------
    # get point estimate to start VI
    #-------------------------------------------------------

    μ₀ = let

        verbose ? @printf("Initialising variational inference.\n") : nothing

        local opt = Optim.Options(show_trace = true, iterations = 3_000, show_every = 25)

        local aux(x) = -helper(f(x))

        optimize(aux, 1*randn(rg, 1*L), NelderMead(), opt).minimizer
     
    end

    
    #-------------------------------------------------------
    # initial covariance root is spherical
    #-------------------------------------------------------
    
    Cdiag = 0.1 * ones(1*L)
    
    
    #-------------------------------------------------------
    # Capture progress
    #-------------------------------------------------------

    tracknegelbo = TrackNegativeElbo(elbo, Inf, zeros(1*L))

    counter = 0

    function callback(_)
        
        counter += 1
        
        if mod(counter, 200) == 1
        
            @printf("Iter %d\t elbo is %f,\t test elbo is %f\n",counter, tracknegelbo.bestsofar_f, -testelbo(tracknegelbo.bestsofar_x))
        
        else
            
            mod(counter, 10) ==1 ? @printf("Iter %d\t elbo is %f\n",counter,tracknegelbo.bestsofar_f) : nothing

        end

        return false

    end
    #-------------------------------------------------------
    # optimise elbo and get optimal variational parameters
    #-------------------------------------------------------

    opt = Optim.Options(show_trace = false, iterations = iterations, show_every=1, callback = callback)
   
    θ = optimize(tracknegelbo, [μ₀; Cdiag], NelderMead(), opt).minimizer
    

    #----------------------------------------------
    # instantiate approximate posterior
    #----------------------------------------------

    q = transformedgaussianposterior(elbo, θ)


    #----------------------------------------------
    # Draw samples from predictive distribution
    #----------------------------------------------
    
    # samplepredict(ttest0::Array{T, 1}) where T<:Real = samplepredict([ttest0 for _ in 1:L])


    # function samplepredict(ttest) 

    #     local θ = rand(q)

    #     local α, b, τ = θ[1:L], θ[L+1:2L], [0; θ[2L+1:3L-1]]


    #     local K = delayedCovariance(kernel, α, τ, ρfixed, tarray)

    #     local KSobsB = K + Sobs

    #     # matrix for replicating elements

    #     local Q✴ = Qmatrix(length.(ttest))

    #     # N × Ntest
        
    #     local kB✴ = delayedCovariance(kernel, α, τ, ρfixed, tarray, ttest)

    #     # Ntest × Ntest
        
    #     local cB = delayedCovariance(kernel, α, τ, ρfixed, ttest)

    #     # full predictive covariance

    #     local Σpred = Symmetric(cB - kB✴' * (KSobsB \ kB✴)) + 1e-8*I

    #     # predictive mean

    #     local μpred = kB✴' * (KSobsB \ (Y - (Q*b))) + (Q✴ * b)

    #     # draw sample and then split it per filter

    #     local Ntest = length.(ttest)

    #     local sample = rand(MvNormal(μpred, Σpred))

    #     return [sample[(sum(Ntest[1:(i-1)])+1):sum(Ntest[1:i])] for i in 1:L]

    # end
    
    # return q, samplepredict
    return q
end
