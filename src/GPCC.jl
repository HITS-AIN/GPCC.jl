module GPCC

    using PyPlot, BlockArrays, Random

    using Optim, Distributions, LinearAlgebra

    using Printf, MiscUtil
    
    using StatsFuns

    using ELBOfy#, ProgressMeter


    include("delayedCovariance.jl")

    include("simulatedata.jl")

    include("util.jl")

    include("gpccfixdelay_marginaliseb.jl")

    include("getprobabilities.jl")

    include("uniformpriordelay.jl")

    include("gpccvi.jl") 

    include("intersection.jl"); export noisyintersection

    include("infercommonlengthscale.jl"); export infercommonlengthscale


    export simulatetwolightcurves, simulatethreelightcurves, simulatefourlightcurves,
           gpcc, getprobabilities, uniformpriordelay, gpccvi


end
