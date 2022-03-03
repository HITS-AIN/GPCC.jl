module GPCC

    using PyPlot, Crayons, BlockArrays, Random

    using Optim, Distributions, LinearAlgebra, MiscUtil

    using Printf, MiscUtil

    using MLBase, StatsFuns

    using ApproximateVI


    include("delayedCovariance.jl")

    include("gpccvi.jl")

    include("simulatedata.jl")

    include("util.jl")

    include("gpcc2.jl")

    export gpcc2vi, simulatedata, gpcc2

end
