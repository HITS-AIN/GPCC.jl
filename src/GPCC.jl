module GPCC

    using PyPlot, Crayons, BlockArrays, Random

    using Optim, Distributions, LinearAlgebra, MiscUtil

    using Printf

    using MLBase, StatsFuns

    using ApproximateVI


    include("delayedCovariance.jl")

    include("gpcc.jl")

    include("Experiments/Synthetic/simulatedata.jl")

    include("util.jl")

    include("gpccvi.jl")

    export gpcc, simulatedata, gpccvi

end
