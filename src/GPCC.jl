module GPCC

    using PyPlot, Random

    using Optim, Distributions, LinearAlgebra, ForwardDiff

    using Printf, MiscUtil
    
    using StatsFuns, Distances

    using ProgressMeter


    include("setup_problem.jl"); export setup_problem
    
    include("newcov.jl")

    include("simulatedata.jl")

    include("util.jl")

    include("gpccfixdelay_marginaliseb.jl"); 
    include("gpccfixdelay_marginaliseb_sparse.jl"); include("inducingpoints.jl")

    include("getprobabilities.jl")

    include("uniformpriordelay.jl")

    # include("gpccvi_b.jl") 

    include("intersection.jl"); export noisyintersection

    include("infercommonlengthscale.jl"); export infercommonlengthscale


    export simulatetwolightcurves, simulatethreelightcurves, simulatefourlightcurves,
    gpcc,gpcc_sparse, getprobabilities, uniformpriordelay#, gpccvi


end
