
<h1 align="center">GPCC.jl</h1>
<p align="center">
  <img width="253" height="165" src=logo.png>
</p>



## ℹ What is this?

It is a julia implementation of the Gaussian Process Cross Correlation (GPCC) method introduced in 

*A Gaussian process cross-correlation approach to time delay estimation for reverberation mapping of active galactic nuclei*.

## 💾 Installation

Apart from cloning, an easy way of using the package is the following:

1 - Add the registry [AINJuliaRegistry](https://github.com/HITS-AIN/AINJuliaRegistry)

2 - Switch into "package mode" with `]` and add the package with
```
add GPCC
```

The package exposes the following functions of interest to the user: 
- `gpcc`
- `simulatetwolightcurves` and `simulatethreelightcurves`,
- `getprobabilities`
- `uniformpriordelay`.

These functions can be queried in help mode in the Julia REPL. 


If installing `GPCC.jl` in an existing Julia environment, there is a chance one may run into dependency problems that prevent installation. In this case, it is advisable to work in a new environment. That is

```
mkdir("myGPCC")
cd("myGPCC")
# press `]` to enter package mode:
(@v1.7) pkg> activate .
```
and use this environment for installing and working with the package.
When restarting Julia, one can re-enter this environment by simply starting Julia in the respective folder ("myGPCC") and using `activate .` in package mode.



## ▶ How to simulate data

Method `simulatedata` can be used to simulate data in 2 arbitrary (non-physical) bands:
```
using GPCC
tobs, yobs, σobs, truedelays = simulatetwolightcurves() # output omitted
```

<p align="center">
  <img src=simulateddata.png>
</p>

A figure like the one above should show up displaying simulated light curves.

It is important to note how the simulated data are organised because function `gpcc` expects the data passed to it to be organised in the exact same way.
First of all, we note that all three returned outputs are vectors whose elements are vectors (i.e. arrays of arrays) and  that they share the same size:
```
typeof(tobs), typeof(yobs), typeof(σobs) 
size(tobs), size(yobs), size(σobs)
```
Each output contains data for 2 bands.
`tobs` contains the observed times. `tobs[1]` contains the observed times for the 1st band, `tobs[2]` for the 2nd band.
Similarly `yobs[1]` contains the flux measurements for the 1st band and `σobs[1]` the error measurements for the 1st band and so on.
We can plot the data pertaining to the 2nd band as an example:

```
using PyPlot # must be indepedently installed, other plotting packages can be used instead
figure()
errorbar(tobs[2], yobs[2], yerr=σobs[2], fmt="o", label="2nd band")
```



## ▶ How to fit a dataset with `gpcc`

Having generated the simulated data, we will now model them with the GPCC model. To that end we use the function `gpcc`. Options for `gpcc` can be queried in help mode.

```
using GPCC

tobs, yobs, σobs, truedelays = simulatetwolightcurves();

# We choose the Matern32 kernel. Other choices are GPCC.OU, GPCC.rbf, GPCC.matern32, GPCC.matern52
# We fit the model for the given the true delays 
# Note that without loss of generality we can always set the delay of the 1st band equal to zero
# The optimisation of the GP hyperparameters runs for a maximum of 1000 iterations.

loglikel, pred, (α, postb, ρ) = gpcc(tobs, yobs, σobs; kernel = GPCC.matern32, delays = truedelays, iterations = 1000, rhomax = 300)
```
The call returns three outputs:
- the marginal log likelihood `loglikel` reached by the optimiser.
- a function `pred` for making predictions.
- a tuple that contains the scaling coefficients $\alpha$, posterior distribution `postb` (of type [MvNormal](https://juliastats.org/Distributions.jl/stable/multivariate/#Distributions.MvNormal)) for shift $b$  and lengthscale $\rho$ of the latent Gaussian process.

We show below how function `pred` can be used both for making predictions and calculating the predictive likelihood.

## ▶ How to make predictions

Having fitted the model to the data, we can now make predictions. We first define the interval over which we want to predict and use `pred`:
```
t_test = collect(0:0.1:20);
μpred, σpred = pred(t_test);
```

Both `μpred` and `σpred` are arrays of arrays. The $l$-th inner array refers to predictions for the $l$-th band, e.g. `μpred[2]` and `σpred[2]` hold respectively the mean prediction and standard deviation of the $2$-band. We plot the predictions for all bands:


<p align="center">
  <img src=simulateddata_predictions.png>
</p>

```
using PyPlot # must be independently installed, other plotting packages can be used instead

colours = ["blue", "orange"] # define colours

for i in 1:2
    plot(t_test, μpred[i], "-", color=colours[i])
    fill_between(t_test, μpred[i] + σpred[i], μpred[i] - σpred[i], color=colours[i], alpha=0.2) # plot uncertainty tube
end

```




## ▶ How to calculate log-likelihood on test data

Suppose we want to calculate the log-likelihood on some new data (test data perhaps):
```
ttest = [[9.0; 10.0; 11.0], [9.0; 10.0; 11.0]]
ytest = [ [6.34, 5.49, 5.38], [13.08, 12.37, 15.69]]
σtest = [[0.34, 0.42, 0.2], [0.87, 0.8, 0.66]]

pred(ttest, ytest, σtest)
```

## ▶ Evaluating a set of candidate delays

Given the simulated data, suppose we would like to evaluate the posterior probability of a set of candidate delays.
Noting that without loss of generality we can always set the delay of the 1st band equal to zero, we define the following grid of delays:
```
candidatedelays = collect(0.0:0.2:20)
```

We use `map` to run `gpcc` on all candidate delays as follows:
```
using GPCC

using PyPlot # we need this for plotting the posterior probabilities, must be independently installed. Other plotting packages can be used instead

tobs, yobs, σobs, truedelays = simulatetwolightcurves();

helper(delay) = gpcc(tobs, yobs, σobs; kernel = GPCC.matern32, delays = [0;delay], iterations = 1000, rhomax = 300)[1] # keep only first output

loglikel = map(helper, candidatedelays)

figure()

plot(candidatedelays, getprobabilities(loglikel))
```

## ▶ Evaluating a set of candidate delays in parallel

One can easily parallelise cross-validation on multiple cores by simply replacing `map` with `pmap`. Before that, one has to make sure that multiple workers are available:
```
using Distributed

addprocs(4) # add four workers. Alternatively start Julia with mulitple workers e.g. julia -p 4

@everywhere using GPCC # make sure GPCC is made available to all workers

@everywhere using ProgressMeter, Suppressor # need to be independently installed

using PyPlot # we need this to plot the posterior probabilities, must be independently installed. Other plotting packages can be used instead

candidatedelays = collect(0.0:0.1:20)

tobs, yobs, σobs, truedelays = simulatetwolightcurves();

# macro @showprogress below reports progress of pmap with a progress bar
# macro @suppress below suppresses terminal messages produced by gpcc

loglikel = @showprogress pmap(candidatedelays) do delay

  @suppress gpcc(tobs, yobs, σobs; kernel = GPCC.matern32, delays = [0;delay], iterations = 1000, rhomax = 300)[1] # keep only first output

end

figure()

plot(candidatedelays, getprobabilities(loglikel))
```


## ▶ Evaluating a set of candidate delays for 3 light curves

We show an example for calculating the posterior for 3 light curves.
Instead of function `simulatetwolightcurves`, we use function `simulatethreelightcurves` to generate 3 synthetic light curves.
We evaluate the delays using a nested `map`:

```
using GPCC

using ProgressMeter, Suppressor # need to be independently installed

using PyPlot # we need this to plot the posterior probabilities, must be independently installed. Other plotting packages can be used instead

candidatedelays = collect(0.5:0.05:6) # use smaller and finer range

tobs, yobs, σobs, truedelays = simulatethreelightcurves();

out = @showprogress map(d2 -> map(d1 -> (@suppress gpcc(tobs, yobs, σobs; kernel = GPCC.matern32, delays = [0;d1;d2], iterations = 1000, rhomax = 300)[1]), candidatedelays), candidatedelays);

posterior = getprobabilities(reduce(vcat, out));

posterior = reshape(posterior, length(candidatedelays), length(candidatedelays));

figure()
subplot(311)
title("joint posterior")
pcolor(candidatedelays,candidatedelays,posterior)
ylabel("lightcurve 2"); xlabel("lightcurve 3")

subplot(312)
title("marginal posterior for lightcurve 2")
plot(candidatedelays,vec(sum(posterior,dims=2)))

subplot(313)
title("marginal posterior for lightcurve 3")
plot(candidatedelays,vec(sum(posterior,dims=1)))
```

We should obtain a joint posterior and marginal posteriors like the ones plotted below:

<p align="center">
  <img src=2Dposterior.png alt="2Dposterior">
</p>

The above computation can be parallelised by starting additional workers and replacing the outer `map` with a `pmap`.



<!---

## ▶ How to decide between candidate delays using `performcv`

Suppose we did not know what the true delays characterising the simulated light curves were.
In this case we would propose a few candidate delays, like e.g.
```
candidatedelays = 0.0:0.1:5.0
```
and subject them to $5$-fold cross-validation as follows:
```
cvresults = map(candidatedelays) do d
  performcv(tobs, yobs, σobs; kernel =  GPCC.matern32, delays = [0;d], iterations = 1000, numberoffolds = 5)
end
```

We obtain approximate posterior probabilities with:
```
post = getprobabilities(cvresults)
figure()
plot(candidatedelays, post, "o-"); xlabel("delays"); ylabel("prob") # PyPlot must be imported
```

<p align="center">
  <img src=delay_vs_prob.png>
</p>


## ▶ How to use `performcv` on multiple cores

One can easily parallelise cross-validation on multiple cores by simply replacing `map` with `pmap`. Before that, one has to make sure that multiple workers are available:
```
using Distributed
addprocs(2) # add two workers. Alternatively start Julia with mulitple workers e.g. julia -p 2
@everywhere using GPCC # make sure GPCC is made available to all workers

cvresults2 = pmap(candidatedelays) do d
  performcv(tobs, yobs, σobs; kernel = GPCC.matern32, delays = [0;d], iterations = 1000, numberoffolds = 5)
end

post2 = getprobabilities(cvresults2)


# Check that the results are approximately the same.
# Note that results will not be exactly identical as the code does not guarantee
# that the same random seeds are used both in parallel and single worker mode
all(post .≈ post2)
```
-->
