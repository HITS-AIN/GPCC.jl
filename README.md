
<h1 align="center">GPCC.jl</h1>
<p align="center">
  <img width="253" height="165" src=logo.png>
</p>





## 💾 Installation

Apart from cloning, an easy way of using the package is the following:

1 - Add the registry [CollaborativeAstronomyJulia](https://github.com/ngiann/CollaborativeAstronomyJulia).

2 - Switch into "package mode" with ```]``` and add the package with
```
add GPCC
```

The package exposes four functions that may be of interest to the user: `gpcc`, `simulatedata`, `getprobabilities` and `performcv`.
These functions can be queried in help mode in the Julia REPL. 

In case you are installing `GPCC.jl` in an existing Julia environment, there is a chance one may run into dependency problems that prevent installation. In this case, it is advisable to work in a new environment. That is

```
mkdir("myGPCC")
cd("myGPCC")
# press `]` to enter package mode:
(@v1.6) pkg> activate .
```
and use this environment for installing and working with the package.
Having exited Julia, one can enter the created environment again by simply starting Julia in the respective folder and using `activate .` in package mode.



## ▶ How to simulate data

Method `simulatedata` can be used to simulate data in 2 arbitrary (non-physical) bands:
```
using GPCC
tobs, yobs, σobs, truedelays = simulatedata() # output omitted
```

<p align="center">
  <img src=simulateddata.png>
</p>

A figure, like the one above, should show up displaying simulated light curves.

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
using PyPlot # must be indepedently installed
errorbar(tobs[2], yobs[2], yerr=σobs[2], fmt="o", label="2nd band")
```



## ▶ How to fit a dataset with `gpcc`

Having generated the simulated data, we will now fit them with the GPCC model. To that end we use the function `gpcc`. Options for `gpcc` can be queried in help mode.

```
using GPCC

tobs, yobs, σobs, truedelays = simulatedata();

# We choose the rbf kernel. Other choices are GPCC.OU / GPCC.rbf / GPCC.matern32.
# We fit the model for the given the true delays 
# Note that without loss of generality we can always set the delay of the 1st band equal to zero.
# The optimisation of the GP hyperparameters runs for a maximum of 1000 iterations.

minopt, pred, posterioroffsetb = gpcc(tobs, yobs, σobs; kernel = GPCC.rbf, delays = truedelays, iterations = 1000)
```
The call returns three outputs:
- the (local) optimum marginal likelihood `minopt` reached by the optimiser.
- a function `pred` for making predictions.
- the posterior distribution of the offset vector `posterioroffsetb` as an object of type [MvNormal](https://juliastats.org/Distributions.jl/stable/multivariate/#Distributions.MvNormal).


## ▶ How to make predictions

Having fitted the model to the data, we can now make predictions. We define the interval over which we want to predict and use `pred`:
```
t_test = collect(0:0.1:20);
μpred, σpred = pred(t_test);
```

Both `μpred` and `σpred` are arrays of arrays. The $l$-th inner array refers to predictions for the $l$-th band, e.g. `μpred[2]` and `σpred[2]` hold respectively the mean prediction and standard deviation of the $2$-band. We plot the predictions for all bands:


<p align="center">
  <img src=simulateddata_predictions.png>
</p>

```
colours = ["blue", "orange"] # define colours

for i in 1:2
    plot(t_test, μpred[i], "-", color=colours[i])
    fill_between(t_test, μpred[i] + σpred[i], μpred[i] - σpred[i], color=colours[i], alpha=0.2) # plot uncertainty tube
end

```




## ▶ How to calculate log-likelihood on test data



## ▶ How to decide between candidate delays using `performcv`

Suppose we did not know what the true delays characterising the simulated light curves were.
In this case we would propose a few candidate delays, like 
```
candidatedelays = 1.0:0.1:3.0
```
and subject them to $5$-fold cross-validation as follows:
```
cvresults = map(candidatedelays) do d
  performcv(tobs, yobs, σobs; kernel = GPCC.rbf, delays = [0;d], iterations = 1000, numberoffolds = 5)
end
```

We obtain approximate posterior probabilities with:
```
post = getprobabilities(cvresults)
plot(candidatedelays, post, "o-"); xlabel("delays"); ylabel("prob") # PyPlot must be imported
```

<p align="center">
  <img src=delay_vs_prob.png>
</p>


## ▶ How to use `performcv` on multiple cores

One can easily parallelise cross-validation on multiple cores by simply replacing `map` with `pmap`. Before that one has to make sure that multiple workers are available:
```
using Distributed
addprocs(2) # add two workers
@everywhere using GPCC # make sure GPCC is made available to all workers

cvresults2 = pmap(candidatedelays) do d
  performcv(tobs, yobs, σobs; kernel = GPCC.rbf, delays = [0;d], iterations = 1000, numberoffolds = 5)
end

post2 = getprobabilities(cvresults2)


# Check that the results are the same.
# Note that results will not be exactly identical as the code does not guarantee
# that the same random seeds are used both in parallel and single worker mode
all(post .≈ post2)
```

## 🔵 Experimental results (THIS WILL BE MOVED TO THE PAPER RELATED PACKAGE)

See [here](https://github.com/ngiann/GPCCExperiments) for experimental results.


## 🔵 Example using a real dataset (THIS WILL BE MOVED TO THE PAPER RELATED PACKAGE)

We use the package [GPCCData](https://github.com/ngiann/GPCCData.jl) to access real observations.


```
using GPCC, GPCCData

# Following packages need to be independently installed. 
# ProgressMeter provides a progress bar while the user waits and Suppressor surpresses output to the terminal
using ProgressMeter, Suppressor 


# load data
tobs, yobs, σobs, _ = readdataset(source="Mrk6")

# Let's look at how data are organised. All of the three arrays have the same structure. They are all arrays of arrays.
display(type(tobs)), display(type(yobs)), display(type(σobs))

# Each array contains 2 inner arrays, one for each observed band (The number of bands is referred to as L in the paper).
length.(tobs)
length.(yobs)
length.(σobs)

# We define an array of candidate delay vectors. Without loss of generalisation, the delay that corresponds to the first light curve is fixed to 0.
delays = [[0;d] for d in 0:0.2:100]

# Check number of candidate delay vectors
length(delays)

# We want to run cross-validation for all candidate delay vectors.

# Do "warmup" first for Julia
@showprogress map(D -> (@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, iterations=1000, numberofrestarts=3, delays = D, kernel = GPCC.matern32)), delays[1:2])

# Do proper run 
out = @showprogress map(D -> (@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, iterations=2000, numberofrestarts=1, delays = D, kernel = GPCC.matern32)), delays)

# estimate posterior probability
getprobabilities(out)
```

## 🔵 Example using a real dataset on multiple cores (THIS WILL BE MOVED TO THE PAPER RELATED PACKAGE)

The above example can be easily be parallelised, i.e. we can try out the candidate delays in parallel.
We need to start julia with multiple processes e.g. "julia -p 16" starts Julia with 16 workers.
Alternatively, we can create more workers within Julia with:
```
using Distributed
addprocs(16) # put here number of available cores
```

We repeat the script from above with minor changes marker with ⚠.
We discard the lines of code inspecting the size and type of the variables.

```
using GPCCData

@everywhere using GPCC  # ⚠ @everywhere makes packages available to all workers ⚠

# Following packages need to be independently installed. 
# ProgressMeter provides a progress bar while the user waits and Suppressor surpresses output to the terminal

@everywhere using ProgressMeter, Suppressor # ⚠ again we use @everywhere ⚠


# load data
tobs, yobs, σobs, _ = readdataset(source="Mrk6")

# We define an array of candidate delay vectors. Without loss of generalisation, the delay that corresponds to the first light curve is fixed to 0.
delays = [[0;d] for d in 0:0.2:100]

# We want to run cross-validation for all candidate delay vectors in parallel!

# Do "warmup" first for Julia
@showprogress pmap(D -> (@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, iterations=1000, numberofrestarts=3, delays = D, kernel = GPCC.matern32)), delays[1:2*nworkers()]) # ⚠ use pmap instead map ⚠

# Do proper run 
out = @showprogress pmap(D -> (@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, iterations=2000, numberofrestarts=1, delays = D, kernel = GPCC.matern32)), delays) # ⚠ use pmap instead map ⚠

# estimate posterior probability
getprobabilities(out)
```
