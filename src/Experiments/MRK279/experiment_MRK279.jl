using ThreadTools

function experiment_MRK279(mass=1e8; iterations=100)

  fixedMdot = 0.01

  @printf("Given mass is %e\n", mass)
  @printf("Given accretion rate is %.2f\n", fixedMdot)


  lambda, tobs, yobs, σobs = readMRK279()


  delays = [shakuradelay(l; M=mass, R=fixedMdot) for l in lambda]

  garray = [x -> boxcar(x; width = w*2) for w in delays]

  callme(_dummy_) = traingp2(tobs, yobs, garray, iterations = iterations)

  results = tmap(callme , 1:5)

  bestindex = argmin([r[1] for r in results])

  fmin, pred = results[bestindex][1], results[bestindex][2]


  xtest = LinRange(minimum(tobs[1]), maximum(tobs[1]), 1000)

  for i in 1:length(yobs)

    plot(xtest, pred(xtest)[1][i],"--")

  end

  return fmin

end
