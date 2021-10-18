# nice ~/julia-1.6.1/bin/julia -O3 -t7
using StatsFuns
using MLBase

include("init.jl")

include("performcv.jl")



include("traingpwithobservednoisefast2.jl")

tobs, yobs, σobs, garray, delays, lambda = simulatedatafromgp(σ=0.5)

# warmup
performcv(tarray=tobs, yarray=yobs, stdarray = σobs, delayarray = [shakuradelay(l; M=1e8, R=0.5) for l in lambda], iterations=2, minimumSqLen = 0.1, numberoffolds=2)

close("all")

mrange = collect(logrange(1e7, 1e8, 20))

results = Vector{Any}(undef, length(mrange))

for index in 1:length(mrange)

    @printf("\n\n %d \t mass=%e\n\n", index, mrange[index])

    results[index] = performcv(tarray=tobs, yarray=yobs, stdarray = σobs, delayarray = [shakuradelay(l; M=mrange[index], R=0.5) for l in lambda], iterations=3_000, minimumSqLen = 0.1, numberoffolds=5)

end
