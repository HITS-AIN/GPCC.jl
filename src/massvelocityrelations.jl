#---------------------------------------------------

function delaygivenvirialmass(; Mv = Mv, V = V)

    c = 2.99792458e10 # speed of light in cm s-1

    g = 6.67259e-8    # Gravitational constant in cm3 gram-1 s-2

    aux = (V * 1000.0 * 100.0)^2 * 86400.0 * c / g

    delay = (Mv * 1.989e33) / aux # in days

    return delay

end


function velocitygivenvirialmass(; Mv = Mv, τ = τ)

    c = 2.99792458e10 # speed of light in cm s-1

    g = 6.67259e-8    # Gravitational constant in cm3 gram-1 s-2

    sqrt(Mv * 1.989e33) / sqrt(τ * 86400.0 * c/g) / (1000 * 100)

end
