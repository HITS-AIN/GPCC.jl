function delayedCovariance(kernel, scale, delays, ρ, x, y)

    @assert(all(scale .> 0))

    if ρ <= 0
        error(@sprintf("ρ=%.8f is <= 0", ρ))
    end

    # number of bands

    L = length(scale) ; @assert(L == length(x) == length(y))

    # number of observations per band

    Nx = map(length, x)

    Ny = map(length, y)

    # it is convenient to use a block array to calculate the covariance matrix

    B = PseudoBlockArray(zeros(eltype(scale), sum(Nx)*L, sum(Ny)*L), Nx, Ny)

    for i in 1:L

        for j in 1:L

            B[Block(i,j)] = [scale[i] * scale[j] * kernel(x₁-delays[i],x₂-delays[j]; ρ=ρ)  for x₁ in x[i], x₂ in y[j]]

        end

    end

    return Matrix(B)

end


delayedCovariance(kernel, scale, delays, ρ, x) = delayedCovariance(kernel, scale, delays, ρ, x, x)
