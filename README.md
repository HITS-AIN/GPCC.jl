
<h1 align="center">GPCC.jl</h1>
<p align="center">
  <img width="253" height="165" src=logo.png>
</p>





## 💾 Installation

Apart from cloning, an easy way of using the package is the following:

1 - Add the registry [CollaborativeAstronomyJulia](https://github.com/ngiann/CollaborativeAstronomyJulia). (❗ change this to [AINJuliaRegistry](https://github.com/HITS-AIN/AINJuliaRegistry) once public)

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
figure()
errorbar(tobs[2], yobs[2], yerr=σobs[2], fmt="o", label="2nd band")
```



## ▶ How to fit a dataset with `gpcc`

Having generated the simulated data, we will now fit them with the GPCC model. To that end we use the function `gpcc`. Options for `gpcc` can be queried in help mode.

```
using GPCC

tobs, yobs, σobs, truedelays = simulatedata();

# We choose the Matern32 kernel. Other choices are OU(), RBF(), Matern32(), Matern52().
# We fit the model for the given the true delays 
# Note that without loss of generality we can always set the delay of the 1st band equal to zero.
# The optimisation of the GP hyperparameters runs for a maximum of 1000 iterations.

minopt, pred, posterioroffsetb, scalingcoeff, lengthscale = gpcc(tobs, yobs, σobs; kernel = Matern32(), delays = truedelays, iterations = 1000)
```
The call returns five outputs:
- the (local) optimum marginal likelihood `minopt` reached by the optimiser.
- a function `pred` for making predictions.
- the posterior distribution of the offset vector `posterioroffsetb` as an object of type [MvNormal](https://juliastats.org/Distributions.jl/stable/multivariate/#Distributions.MvNormal).
- scaling coefficients
- lengthscale of latent Gaussian process

We show below that function `pred` can be used both for making predictions and calculating the predictive likelihood.

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

Suppose we want to calculate the log-likelihood on some new data (test data perhaps):
```
ttest = [[9.0; 10.0; 11.0], [9.0; 10.0; 11.0]]
ytest = [ [6.34, 5.49, 5.38], [13.08, 12.37, 15.69]]
σtest = [[0.34, 0.42, 0.2], [0.87, 0.8, 0.66]]

pred(ttest, ytest, σtest)
```



## ▶ How to decide between candidate delays using `performcv`

Suppose we did not know what the true delays characterising the simulated light curves were.
In this case we would propose a few candidate delays, like 
```
candidatedelays = 0.0:0.1:5.0
```
and subject them to $5$-fold cross-validation as follows:
```
cvresults = map(candidatedelays) do d
  performcv(tobs, yobs, σobs; kernel = Matern32(), delays = [0;d], iterations = 1000, numberoffolds = 5)
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
  performcv(tobs, yobs, σobs; kernel = Matern32(), delays = [0;d], iterations = 1000, numberoffolds = 5)
end

post2 = getprobabilities(cvresults2)


# Check that the results are approximately the same.
# Note that results will not be exactly identical as the code does not guarantee
# that the same random seeds are used both in parallel and single worker mode
all(post .≈ post2)
```
