function delayedCovariance(scale, delays, ℓ², x, y)

    @assert(all(scale .> 0))

    @assert(ℓ² > 0)

    # number of bands

    L = length(scale) ; @assert(L == length(x) == length(y))

    # number of observations per band

    Nx = map(length, x)

    Ny = map(length, y)

    # it is convenient to use a block array to calculate the covariance matrix

    B = PseudoBlockArray(zeros(eltype(scale), sum(Nx)*L, sum(Ny)*L), Nx, Ny)

    for i in 1:L

        for j in 1:L

            B[Block(i,j)] = [scale[i] * scale[j] * rbf(x₁-delays[i],x₂-delays[j]; ℓ²=ℓ²)  for x₁ in x[i], x₂ in y[j]]

        end

    end

    return Matrix(B)

end


delayedCovariance(scale, delays, ℓ², x) = delayedCovariance(scale, delays, ℓ², x, x)
