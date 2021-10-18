function experiment2_b(; seed = 1, σ = 1.0)


    tobs, yobs, σobs, g = simulatedatafromgp(seed = seed, σ = σ)


    candidatesg = [[x -> boxcar(x; width = w) for w in [0.1, 1.0, 1.5]],
                   [x -> boxcar(x; width = w) for w in [0.2, 1.1, 1.6]],
                   [x -> boxcar(x; width = w) for w in [0.3, 1.2, 1.7]],
                   [x -> boxcar(x; width = w) for w in [0.1, 1.0, 1.2]]]


    map( g -> traingpwithobservednoise(tobs, yobs, σobs, g)[1], candidatesg)


end
