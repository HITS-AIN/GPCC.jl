function setup_problem(tobs, yobs, σobs; kernel = kernel, iterations = 1_000, numberofrestarts = 10, initialrandom = 10, ρmin = 0.1, ρmax = 300.0, rng = AbstractRNG=Random.GLOBAL_RNG)

    L = length(tobs)

    ρfixed = infercommonlengthscale(tobs, yobs, σobs; kernel = kernel, iterations = iterations, rng = rng, numberofrestarts = numberofrestarts, initialrandom = initialrandom, ρmin = ρmin, ρmax = ρmax, verbose = false)[3][3]

    function helper(delay::T...) where T<:Real
  
        gpcc(tobs, yobs, σobs; rng = rng, kernel = kernel, delays = vcat(0, delay...), iterations = iterations, rhomax = ρmax, ρfixed = ρfixed)[1]
    
    end

    function helper(delay::Vector{T}...) where T<:Real
  
        @showprogress tmap1(x->helper(x...), Iterators.product(delay...));

    end


    function helper(delays::Vector{T}) where T<:Real
  
        if L == 2

            @showprogress tmap1(x -> helper(x), delays)

        else

            helper(ntuple(i -> delays, L-1)...)

        end
    end

    function helper(delays::AbstractRange{T}) where T<:Real
    
        @showprogress tmap1(x -> helper(x), delays)

    end

end