function experiment1(; seed = 1, σ = 1e-1)

    # generate mock data

    tobs, yobs, σobs, g = simulatedatafromfftconvolution(seed = seed, σ = σ)

    # fit GP

    fmin, pred = traingpwithobservednoise(tobs, yobs, σobs, g)

    # make predictions

    xtest = LinRange(minimum(reduce(vcat, tobs)), maximum(reduce(vcat, tobs)), 500)

    μ, σ = pred(xtest)

    # plot results

    colours = ["r", "g", "b"]

    figure(-1)

    title("lines of the same colour must overlap")

    for i in 1:3

        plot(xtest, μ[i], colours[i]*"--")

        fill_between(xtest, μ[i] .+ 2σ[i], μ[i] .- 2σ[i], alpha=0.1, color=colours[i])

    end

    nothing

end
