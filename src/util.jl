function informuser(; seed = seed, iterations = iterations, initialrandom = initialrandom, numberofrestarts = numberofrestarts, JITTER = JITTER, ρmin = ρmin, ρmax = ρmax, Σb = Σb)

    colourprint(@sprintf("Running with random seed %d\n", seed), foreground = :light_blue, bold = true)
    @printf("\t iterations             = %d\n", iterations)
    @printf("\t initialrandom          = %d\n", initialrandom)
    @printf("\t numberofrestarts       = %d\n", numberofrestarts)
    @printf("\t JITTER                 = %e\n", JITTER)
    @printf("\t ρmin                   = %f\n", ρmin)
    @printf("\t ρmax                   = %f\n", ρmax)
    @printf("\t Σb                     = "); map(x->@printf("%.3f ",x), diag(Σb)); @printf("\n")
end

#---------------------------------------------------

function OU(xᵢ,xⱼ ; ρ=1.0)

    ℓ = sqrt(ρ)

    r = norm(xᵢ - xⱼ)

    exp(-r/ℓ)

end


#---------------------------------------------------

rbf(xᵢ,xⱼ ; ρ=1.0) = exp(-0.5*(xᵢ-xⱼ)^2/(2ρ))

#---------------------------------------------------

function matern32(xᵢ,xⱼ ; ρ=1.0)

    ℓ = sqrt(ρ)

    r = norm(xᵢ - xⱼ)

    (1 + sqrt(3)*r/ℓ) * exp(- sqrt(3)*r/ℓ)

end

#---------------------------------------------------

function matern52(xᵢ,xⱼ ; ρ=1.0)

    ℓ = sqrt(ρ)

    r = norm(xᵢ - xⱼ)

    (1 + sqrt(5)*r/ρ +( 5*r^2)/(3*ℓ^2)) * exp(- sqrt(5)*r/ℓ)

end

#---------------------------------------------------

function Qmatrix(Narray)

    L = length(Narray) # number of time series

    local Q = zeros(sum(Narray), L)

    for l in 1:L

        Q[1+sum(Narray[1:l-1]):sum(Narray[1:l]) ,l] .= 1.0

    end

    Q

end

#---------------------------------------------------

function Qvector(Narray, entries)

    L = length(Narray) # number of time series

    @assert(L == length(entries))

    local Q = zeros(sum(Narray))

    for l in 1:L

        Q[1+sum(Narray[1:l-1]):sum(Narray[1:l])] .= entries[l]

    end

    Q

end

#---------------------------------------------------

function internalversion()
    @printf("newaltbranch 1 - new init mechanism")
end
