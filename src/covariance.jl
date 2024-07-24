function delayedCovariance(kernel, scale, delays, ρ, x, y)

    ρ <= 0 ? error(@sprintf("ρ=%.8f is <= 0", ρ)) : nothing

    # number of bands

    @assert(length(scale) == length(x) == length(y))

    delayedx = reduce(vcat, [x.-d for (x,d) in zip(x, delays)])

    delayedy = reduce(vcat, [y.-d for (y,d) in zip(y, delays)])

    Sx = Diagonal(Qmatrix(length.(x))*scale)
    
    Sy = Diagonal(Qmatrix(length.(y))*scale)

    Sx*[kernel(x,y;ρ=ρ) for x in delayedx, y in delayedy]*Sy

end


function delayedCovariance(kernel, scale, delays, ρ, x) 
    
    ρ <= 0 ? error(@sprintf("ρ=%.8f is <= 0", ρ)) : nothing

    # number of bands

    @assert(length(scale) == length(x))

    delayedx = reduce(vcat, [x.-d for (x,d) in zip(x, delays)])

    Sx = Diagonal(Qmatrix(length.(x))*scale)

    Sx*[kernel(x,y;ρ=ρ) for x in delayedx, y in delayedx]*Sx

end


covariance_unit_amplitude(kernel, ρ, X)  = Symmetric(covariance_unit_amplitude(kernel, ρ, X, X))

function covariance_unit_amplitude(kernel, ρ, X, Y) 

    ρ <= 0 ? error(@sprintf("ρ=%.8f is <= 0", ρ)) : nothing

    [kernel(x, y; ρ = ρ) for x in X, y in Y]

end