function posteriordelay(tobs, yobs, σobs, candidatedelays; JITTER = 1e-10, kernel = kernel, iterations = 1_000, numberofrestarts = 10, initialrandom = 10, ρmin = 0.1, ρmax = 300.0, rng = AbstractRNG=Random.GLOBAL_RNG)

    L = length(tobs)
    
    ρfixed = infercommonlengthscale(tobs, yobs, σobs, ; kernel = kernel, iterations = iterations, rng = rng, numberofrestarts = numberofrestarts, initialrandom = initialrandom, ρmin = ρmin, ρmax = ρmax, verbose = false)
 
    @printf("Inferred length scale is ρ=%f\n", ρfixed)


    function helper(idx, delay)

        gpcc(tobs[idx], yobs[idx], σobs[idx]; rng = rng, kernel = kernel, delays = vcat(0, delay), iterations = iterations, ρfixed = ρfixed, JITTER = JITTER)[1]
        
    end


    # Define functions for pairs of lightcurves

    functionpairs = [@memoize ThreadSafeDict x -> helper([i;j], x) for i in 1:L for j in (i+1):L]

    
    function H(x...)

        @assert(length(x) == L-1)

        local aux = 0.0

        for i in 1:L-1
            
            aux += functionpairs[i](x[i])

        end

        local counter = L-1

        for i in 2:L, j in (i+1):L
                
            counter += 1
                
            aux += functionpairs[counter](x[j-1]-x[i-1])

        end

        @assert(counter == length(functionpairs))

        return aux

    end

    loglikel = @showprogress tmap1(x -> H(x...), Iterators.product(ntuple(i -> candidatedelays, L-1)...))

    getprobabilities(loglikel)

end