[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub](https://img.shields.io/github/license/HITS-AIN/GPCC.jl)
<a href="https://ascl.net/2303.006"><img src="https://img.shields.io/badge/ascl-2303.006-blue.svg?colorB=262255" alt="ascl:2303.006" /></a>


<h1 align="center">GPCC.jl</h1>
<p align="center">
  <img width="253" height="165" src=logo.png>
</p>

The repo is currently being updated. Please be aware that certain functionality may not be available in the meanwhile.

## ℹ What is this?

This is a Julia implementation of the Gaussian Process Cross Correlation (GPCC) method introduced in 

[*A Gaussian process cross-correlation approach to time delay estimation for reverberation mapping of active galactic nuclei*](https://github.com/HITS-AIN/GPCCpaper).

GPCC is a probabilistic alternative to the Interpolated Cross Correlation function (ICCF). Advantages over the ICCF include:

- Outputs a probability distribution for the delay.
- It can incorporate a prior on the delay.
- Delivers predictions for out-of-sample data.

## 💾 Installation

Apart from cloning, an easy way of using the package is the following:

1 - Add the registry [AINJuliaRegistry](https://github.com/HITS-AIN/AINJuliaRegistry)

2 - Switch into "package mode" with `]` and add the package with
```
add GPCC
```

The package exposes the following functions of interest to the user: 
- `simulatetwolightcurves` and `simulatethreelightcurves`,
- `infercommonlengthscale`,
- `gpcc`,
- `uniformpriordelay`.

These functions can be queried in help mode in the Julia REPL. 


## 🚀 An important note about performance

*(This note is not specific to the GPCC package; it applies in general whenever BLAS threads run concurrently to julia threads.)*

The package supports the parallel evaluation of candidate delays.
To that end, start julia with multiple threads. For instance, you can start julia with 8 threads using `julia -t8`.
We recommend to use as many threads as physical cores.

To get the most performance, please read this note [here](https://carstenbauer.github.io/ThreadPinning.jl/dev/explanations/blas/) concerning issues when running multithreaded code that makes use of BLAS calls. In most cases, the following instructions suffice:
```
using LinearAlgebra
BLAS.set_num_threads(1) # Set the number of BLAS threads to 1

using ThreadPinning # must be indepedently installed
pinthreads(:cores) # allows you to pin Julia threads to specific CPU-threads 
```

Unless you are using the Intel MKL, we recommend to always use the above code before estimating delays.


## ▶ How to simulate data

Method `simulatetwolightcurves` can be used to simulate data in 2 arbitrary (non-physical) bands:
```
using Plots # must be independently installed,alternative plotting packages may be used. 
using GPCC
tobs, yobs, σobs, truedelays = simulatetwolightcurves() # output omitted

scatter(tobs[1],yobs[1],grid=false,yerror=σobs[1], color="blue")
scatter!(tobs[2],yobs[2],grid=false,yerror=σobs[2], color="orange")
```

<p align="center">
  <img src=simulateddata.png>
</p>

A figure like the one above should show up displaying simulated light curves.

It is important to note how the simulated data are organised because function `gpcc` expects the data passed to it to be organised in the exact same way.
First of all, we note that all three returned outputs are vectors whose elements are vectors (i.e. arrays of arrays) and  that they share the same size:
```
# try this out in repl
typeof(tobs), typeof(yobs), typeof(σobs) 
size(tobs), size(yobs), size(σobs)
```
Each output contains data for 2 bands.
`tobs` contains the observed times. `tobs[1]` contains the observed times for the 1st band, `tobs[2]` for the 2nd band.
Similarly `yobs[1]` contains the flux measurements for the 1st band and `σobs[1]` the error measurements for the 1st band and so on.


## ▶ How to fit a dataset with `gpcc`

Having generated the simulated data, we will now model them with the GPCC model. To that end we use the function `gpcc`. Options for `gpcc` can be queried in help mode.

```
using GPCC

tobs, yobs, σobs, truedelays = simulatetwolightcurves();

# We first determine the lengthscale for the GPCC with the following call.
# We choose the rbf kernel. Other choices are GPCC.OU, GPCC.matern32, GPCC.matern52

ρ = infercommonlengthscale(tobs, yobs, σobs; kernel = GPCC.rbf, iterations = 1000)

# We fit the model for the given the true delays and the estimate lengthscale ρ we got in the line above.
# Note that without loss of generality we can always set the delay of the 1st band equal to zero
# The optimisation of the GP hyperparameters runs for a maximum of 1000 iterations.

loglikel, α, postb, pred = gpcc(tobs, yobs, σobs; kernel = GPCC.rbf, delays = truedelays, iterations = 1000, ρfixed = ρ)
```
The call returns three outputs:
- the marginal log likelihood `loglikel` reached by the optimiser,
-  the scaling coefficients $\alpha$,
-  posterior distribution `postb` (of type [MvNormal](https://juliastats.org/Distributions.jl/stable/multivariate/#Distributions.MvNormal)) for shift $b$,
- a function `pred` for making predictions.

We show below how function `pred` can be used both for making predictions and calculating the predictive likelihood.

## ▶ How to make predictions

Having fitted the model to the data, we can now make predictions. We first define the interval over which we want to predict and use `pred`:
```
t_test = collect(-3:0.1:65);
μpred, σpred = pred(t_test);
```

Both `μpred` and `σpred` are arrays of arrays. The $l$-th inner array refers to predictions for the $l$-th band, e.g. `μpred[2]` and `σpred[2]` hold respectively the mean prediction and standard deviation of the $2$-band. We plot the predictions for all bands:


<p align="center">
  <img src=simulateddata_predictions.png>
</p>

```
using Plots # must be independently installed, other plotting packages can be used instead

plot(t_test, μpred[1]; ribbon = σpred[1], fillalpha=0.2, color="blue")
plot!(t_test, μpred[2]; ribbon = σpred[2], fillalpha=0.2, color="orange")
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

using Plots # we need this for plotting the posterior probabilities, must be independently installed. Other plotting packages can be used instead

tobs, yobs, σobs, truedelays = simulatetwolightcurves();

# Get estimate for lengthscale parameter ρ
ρ = infercommonlengthscale(tobs, yobs, σobs; kernel = GPCC.rbf, iterations = 1000)

helper(delay) = gpcc(tobs, yobs, σobs; kernel = GPCC.rbf, delays = [0;delay], iterations = 1000, ρfixed = ρ)[1] # keep only first output

loglikel = map(helper, candidatedelays)

plot(candidatedelays, getprobabilities(loglikel))
```

## ▶ Evaluating a set of candidate delays in parallel threads

One can easily parallelise the posterior estimation by simply replacing `map` with `tmap` provided by the package [ThreadTools.jl](https://github.com/baggepinnen/ThreadTools.jl). 
Package `ThreadTools.jl` needs to be independently installed. Before that, one has to make sure that multiple threads are available by starting Julia with e.g. `julia -t 4` option:
```
using Plots # we need this to plot the posterior probabilities, must be independently installed. Other plotting packages can be used instead

using GPCC

using ProgressMeter, ThreadTools # need to be independently installed

candidatedelays = collect(0.0:0.1:20);

tobs, yobs, σobs, truedelays = simulatetwolightcurves();

ProgressMeter.ncalls(::typeof(tmap), ::Function, args...) = ProgressMeter.ncalls_map(args...) # this line makes ProgressBar work with ThreadTools, see https://github.com/timholy/ProgressMeter.jl#adding-support-for-more-map-like-functions


# Get estimate for lengthscale parameter ρ
ρ = infercommonlengthscale(tobs, yobs, σobs; kernel = GPCC.rbf, iterations = 1000)

helper(delay) = gpcc(tobs, yobs, σobs; kernel = GPCC.rbf, delays = [0;delay], iterations = 1000, ρfixed = ρ)[1] # keep only first output

loglikel = @showprogress tmap(helper, candidatedelays)

plot(candidatedelays, getprobabilities(loglikel))
```


## ▶ Evaluating a set of candidate delays for 3 light curves in parallel

We show an example for calculating the posterior for 3 light curves in parallel.
To do this, we need to start Julia with multiple threads, e.g. `julia -t 4` starts Julia with 4 threads.
Instead of function `simulatetwolightcurves`, we use function `simulatethreelightcurves` to generate 3 synthetic light curves.
We evaluate the delays using a `tmap` inside a `map`:

```
using Plots # we need this to plot the posterior probabilities, must be independently installed. Other plotting packages can be used instead

using GPCC

using ProgressMeter, ThreadTools # need to be independently installed

ProgressMeter.ncalls(::typeof(tmap), ::Function, args...) = ProgressMeter.ncalls_map(args...)
# this line makes ProgressBar work with ThreadTools, see https://github.com/timholy/ProgressMeter.jl#adding-support-for-more-map-like-functions

candidatedelays = collect(0.5:0.05:6) # use smaller and finer range

tobs, yobs, σobs, truedelays = simulatethreelightcurves();

# Get estimate for lengthscale parameter ρ
ρ = infercommonlengthscale(tobs, yobs, σobs; kernel = GPCC.rbf, iterations = 1000)

out = @showprogress map(d2 -> tmap(d1 -> (gpcc(tobs, yobs, σobs; kernel = GPCC.rbf, delays = [0;d1;d2], iterations = 1000, ρfixed = ρ)[1]), candidatedelays), candidatedelays);

posterior = getprobabilities(reduce(vcat, out));
posterior = reshape(posterior, length(candidatedelays), length(candidatedelays));

p1 = contour(candidatedelays, candidatedelays, posterior, title = "joint posterior", ylabel="lightcurve 2", xlabel="lightcurve 3", aspect_ratio = :equal, fontsize=10,  colorbar=false)
p2 = plot(candidatedelays, vec(sum(posterior,dims=2)), title="marginal for lightcurve 2", label = false, fontsize=10)
p3 = plot(candidatedelays, vec(sum(posterior,dims=1)), title="marginal for lightcurve 3", label = false, fontsize=10)

plot(p2,p1,plot(legend=false,grid=false,foreground_color_subplot=:white),p3,layout=(2,2))   
```

We should obtain a joint posterior and marginal posteriors similar to the ones plotted below:

<p align="center">
  <img src=2Dposterior.png alt="2Dposterior">
</p>


❗ Running GPCC on three light curves can be a very lengthy computation! This is because GPCC will try out in a brute force manner all possible delay combinations. Package (FasterGPCC.jl)[https://github.com/HITS-AIN/FasterGPCC.jl] addresses this issue to a certain extend.
