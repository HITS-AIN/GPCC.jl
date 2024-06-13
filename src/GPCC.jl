module GPCC

    using PyPlot, BlockArrays, Random

    using Optim, Distributions, LinearAlgebra, StatsFuns

    using Printf, MiscUtil #, Suppressor

    using StatsFuns # MLBase

    using ELBOfy


    include("delayedCovariance.jl")

    include("simulatedata.jl")

    include("util.jl")


    # include("UNUSED/gpccvi.jl")

    include("gpccfixdelay_marginaliseb.jl")

    # include("performcv.jl")

    include("getprobabilities.jl")

    include("uniformpriordelay.jl")

    include("gpccvi.jl"); export gpccvi


    export simulatetwolightcurves, simulatethreelightcurves, simulatefourlightcurves,
           gpcc, getprobabilities, uniformpriordelay


end
