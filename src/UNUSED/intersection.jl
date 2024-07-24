#######################################################################
function noisyintersection(; posterior = posterior, u = u)
#######################################################################
    
    D = length(u)

    @assert(norm(u) ≈ 1.0)

    # hack: we care only about the first 2D arguments, the rest don't matter.
    # The correct way to deal with would have been to partition the Gaussian
    # and the apply only the relevant part of the non linear transformation.

    _x = mean([rand(posterior) for _ in 1:30])[2D+1:end]

    logq(x) = logpdf(posterior, [x;_x])


    #----------------------------------------
    function unpack(x)
    #----------------------------------------

        @assert(length(x) == D + 1)

        local α = softplus.(x[1:D])

        local t = x[D+1]

        return α, t

    end

    
    #----------------------------------------
    function logp(x₀)
    #----------------------------------------

        REPEATS = 10

        objective(α, t) = logq([α; u*x₀ - α*t])

        opt = Optim.Options(iterations = 10_000, show_trace = false, show_every = 1)

        res = [optimize(x -> -objective(unpack(x)...), 3*randn(D+1),  NelderMead(), opt) for _ in 1:REPEATS]

        bestindex = argmin([r.minimum for r in res])

        -res[bestindex].minimum # return log likelihood of fitted line

    end


    return logp

end
    