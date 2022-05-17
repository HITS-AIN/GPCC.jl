module GPCC

    using PyPlot, BlockArrays, Random

    using Optim, Distributions, LinearAlgebra, StatsFuns

    using Printf, MiscUtil, Suppressor

    using MLBase, StatsFuns

    # using ApproximateVI


    include("delayedCovariance.jl")

    include("simulatedata.jl")

    include("util.jl")

    include("gpccfixdelay_verifications.jl")

    include("gpccfixdelay.jl")

    include("performcv.jl")

    include("getprobabilities.jl")


    export simulatedata, gpcc, performcv, getprobabilities


end
