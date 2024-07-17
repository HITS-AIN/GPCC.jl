module GPCC

    using PyPlot, Random

    using Optim, Distributions, LinearAlgebra, ForwardDiff

    using Printf, MiscUtil
    
    using StatsFuns, Distances

    # using ELBOfy#, ProgressMeter


    include("newcov.jl")

    include("simulatedata.jl")

    include("util.jl")

    include("gpccfixdelay_marginaliseb.jl")

    include("getprobabilities.jl")

    include("uniformpriordelay.jl")

    # include("gpccvi_b.jl") 

    include("intersection.jl"); export noisyintersection

    include("infercommonlengthscale.jl"); export infercommonlengthscale


    export simulatetwolightcurves, simulatethreelightcurves, simulatefourlightcurves,
           gpcc, getprobabilities, uniformpriordelay#, gpccvi


end
