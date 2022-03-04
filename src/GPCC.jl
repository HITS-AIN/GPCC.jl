module GPCC

    using PyPlot, Crayons, BlockArrays, Random

    using Optim, Distributions, LinearAlgebra, StatsFuns

    using Printf, MiscUtil, Suppressor

    using MLBase, StatsFuns

    using ApproximateVI


    include("delayedCovariance.jl")

    include("gpccvi.jl")

    include("simulatedata.jl")

    include("util.jl")

    include("gpcc2.jl")

    include("gpccfixdelay.jl")

    include("performcv.jl")

    include("getprobabilities.jl")

    include("massvelocityrelations.jl")

    export gpcc2vi, simulatedata, gpcc2, gpccfixdelay, performcv, getprobabilities,
        delaygivenvirialmass, velocitygivenvirialmass

end
