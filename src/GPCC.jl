module GPCC

    using Random # PyPlot, 

    using Optim, Distributions, LinearAlgebra, ForwardDiff

    using Printf, MiscUtil
    
    using StatsFuns, Distances

    using ProgressMeter
    

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
    export simulatetwolightcurves, simulatethreelightcurves
    export gpcc
    export getprobabilities, uniformpriordelay

end
