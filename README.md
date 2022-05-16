# Gaussian Process Cross Correlation

Julia implementation of Gaussian Process Cross Correlation model.


## How to use

Apart from cloning, an easy way of using the package is the following:

1 - Add the registry [CollaborativeAstronomyJulia](https://github.com/ngiann/CollaborativeAstronomyJulia).

2 - Switch into "package mode" with ```]``` and add the package with
```
add GPCC
```

The package exposes four functions that may be of interest to the user: `gpcc`, `simulatedata`, `getprobabilities` and `performcv`.
These functions can be queried in help mode at the Julia REPL. 

In case you are installing ProbabilisticFluxVariationGradient to an existing Julia environment, there is a chance one may run into dependency problems that prevent installation. In this case, it is advisable to work in a new environment. That is

```
mkdir("myGPCC")
cd("myGPCC")
# press `]` to enter package mode:
(@v1.6) pkg> activate .
```
and use this environment for installing and working with the package.
Having exited Julia, one can enter the created environment again by simply starting Julia in the respective folder and using `activate .` in package mode.


## Experimental results

See [here](https://github.com/ngiann/GPCCExperiments) for experimental results.


## Example using a real dataset


```
tobs, yobs, σobs = simulatedata(seed=1, σ=1, N=[1;1;1]*75, ρ=5);

# try true delays
q = performcv(tarray=tobs, yarray=yobs, stdarray=σobs, iterations=1000, numberofrestarts=7, delays = [0;2;6], kernel = GPCC.matern32);

# try perturbed delay
q2 = performcv(tarray=tobs, yarray=yobs, stdarray=σobs, iterations=1000, numberofrestarts=7, delays = [0;2.3;5.8], kernel = GPCC.matern32);

# estimate posterior probability
getprobabilities([q,q2])

# output is:
#2-element Vector{Float64}:
# 0.6356097291597814
# 0.3643902708402206

```
