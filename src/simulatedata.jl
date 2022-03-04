#############################################################################
function simulatedata(; σ = 0.1, seed = 1, N = [50; 40; 30], ρ = 1.75)
#############################################################################

    rg = MersenneTwister(seed)

    #---------------------------------------------------------------------
    # Define GP parameters
    #---------------------------------------------------------------------

    @show delays = [0.0; 2.0; 6.0]

    @show scale  = [1; 2; 0.5]

    @show shift  = [5; 6.0; 9]


    #---------------------------------------------------------------------
    #
    #---------------------------------------------------------------------

    aux  = [rand(rg, N[l]).> 0.5 for l in 1:length(delays)]


    function samplefrominterval(a)

        draw(x) = x==0 ? rand(rg, Uniform(-10.0, 6.0)) : rand(rg, Uniform(8.0, 25.0))

        return [draw(aᵢ) for aᵢ in a]

    end


    t = [samplefrominterval(aux[l]) for l in 1:length(delays)]


    #---------------------------------------------------------------------
    # Define Gaussian process to draw noisy targets
    #---------------------------------------------------------------------

    C = delayedCovariance(matern32, scale, delays, ρ, t)

    let

        U, S, V = svd(C)

        C = U * Diagonal(max.(1e-6, abs.(S))) * U'

        makematrixsymmetric!(C)

    end


    #---------------------------------------------------------------------
    # Draw targets and arrange in array
    #---------------------------------------------------------------------

    Y = rand(rg, MvNormal(zeros(sum(N)), C))

    y = Vector{Vector{Float64}}(undef, length(delays))

    mark = 0

    for i in 1:length(delays)

        y[i] = Y[mark+1:mark+N[i]] * scale[i] .+ shift[i] .+ σ*randn(rg, N[i])

        mark += N[i]

    end

    figure(0) ; cla()

    for i in 1:length(delays)

      plot(t[i], y[i], "o", label = @sprintf("delay = %.3f", delays[i]))

    end

    legend()

    return t, y, σ


end
