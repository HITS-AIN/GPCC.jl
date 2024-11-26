module GPCC

    using Random # PyPlot, 

    using Optim, Distributions, LinearAlgebra, ForwardDiff

    using Printf, MiscUtil
    
    using StatsFuns, Distances

    using ProgressMeter#, ThreadTools

    # using ELBOfy

    # using Memoization, ThreadSafeDicts

    
    # Following lines makes ProgressMeter work with tmap1

    # ProgressMeter.ncalls(::typeof(tmap1), ::Function, args...) = ProgressMeter.ncalls_map(args...)


    # include("posteriordelay.jl")
    
    include("covariance.jl")

    include("simulatedata.jl")

    include("kernels.jl")

    include("Qmatrix.jl")

    include("Qvector.jl")

    include("gpccfixdelay_marginaliseb.jl"); 
    
    include("getprobabilities.jl")

    include("uniformpriordelay.jl")

    # include("gpccvi.jl") 
    
    include("gp_commonlengthscale.jl")
    
    export infercommonlengthscale
    export posteriordelay
    export simulatetwolightcurves, simulatethreelightcurves#, simulatefourlightcurves, simulatefivelightcurves
    export gpcc
    export getprobabilities, uniformpriordelay
    # export gpccvi


end
