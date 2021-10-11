
softplus(x) = log(1.0 + exp(x))

invsoftplus(y) = log(exp(y) - 1.0)

makepositive(x) = softplus(x)

invmakepositive(x) = invsoftplus(x)

#---------------------------------------------------

boxcar(x; width = width) =  x >= 0.0 && x <= width ? 1.0/width : 0.0

#---------------------------------------------------

logrange(x1, x2, n) = (10^y for y in range(log10(x1), log10(x2), length=n))

#---------------------------------------------------

heaviside(t) = 0.5 * (sign(t) + 1)

#---------------------------------------------------

rbf(xᵢ,xⱼ ; ℓ²=1.0) = exp(-0.5*(xᵢ-xⱼ)^2/(2ℓ²))

#---------------------------------------------------

sigmoid(x) = 1.0 / (1.0 + exp(-x))

unsigmoid(y) = -log(1.0 / y - 1.0)

#---------------------------------------------------

function transformbetween(x, a, b)
    @assert(a<b)
    sigmoid(x)*(b-a) + a
end

function untransformbetween(y, a, b)
    @assert(a<b)
    unsigmoid((y-a)/(b-a))
end

#---------------------------------------------------

function delaygivenvirialmass(; Mv = Mv, V = V)

    c = 2.99792458e10 # speed of light in cm s-1

    g = 6.67259e-8    # Gravitational constant in cm3 gram-1 s-2

    aux = (V * 1000.0 * 100.0)^2 * 86400.0 * c / g

    delay = (Mv / 1.989e33) / aux # in days

    return delay

end
