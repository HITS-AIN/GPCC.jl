function setup_problem(tobs, yobs, σobs; kernel = GPCC.matern52, iterations = 1000, seed = 1, numberofrestarts = 10, initialrandom = 10, ρmin = 0.1, ρmax = 300.0)

    ρfixed = infercommonlengthscale(tobs, yobs, σobs; kernel = kernel, iterations = iterations, seed = seed, numberofrestarts = numberofrestarts, initialrandom = initialrandom, ρmin = ρmin, ρmax = ρmax, verbose = false)[3][3]

    function helper(delay::T...) where T<:Real
  
        gpcc(tobs, yobs, σobs; kernel = kernel, delays = vcat(0, delay...), iterations = iterations, rhomax = ρmax, ρfixed = ρfixed)[1]
    
    end

    function helper(delay::Vector{T}...) where T<:Real
  
        @showprogress tmap(x->helper(x...), Iterators.product(delay...));

    end

    L = length(tobs)

    function helper(delay::Vector{T}) where T<:Real
  
        helper(ntuple(i -> delay, L-1)...)

    end

end