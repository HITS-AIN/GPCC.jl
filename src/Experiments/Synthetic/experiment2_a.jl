function experiment2_a(; seed = 1, σ = 1.0, iterations = 1, N = 50)

    # generate mock data

    tobs, yobs, σobs, _, delays, _ = simulatedatafromgp(seed = seed, σ = σ, N = N)

    # fit GP

    fmin, pred = traingpwithobservednoise(tobs, yobs, σobs, delays; initl2 = 0.75, iterations = iterations)

    # make predictions

    xtest = LinRange(minimum(reduce(vcat, tobs)), maximum(reduce(vcat, tobs)), 500)

    μ, σ = pred(xtest)

    # plot data

    figure(-1) ; cla()

    colours = ["r", "g", "b"]

    for i in 1:3

        plot(tobs[i], yobs[i], colours[i]*"o")

    end

    # plot predictions

    for i in 1:3

        plot(xtest, μ[i], colours[i]*"--")

        fill_between(xtest, μ[i] .+ 2σ[i], μ[i] .- 2σ[i], alpha=0.1, color=colours[i])

    end

    nothing

end
