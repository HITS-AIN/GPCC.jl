# Gaussian Process Cross Correlation

A Julia implementation of the Gaussian Process Cross Correlation model.


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

We use the package [GPCCVirialDatasets](https://github.com/ngiann/GPCCVirialDatasets.jl) to access real observations.


```
using GPCC, GPCCVirialDatasets

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

## Example using a real dataset on multiple cores

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
@everywhere using GPCC, GPCCVirialDatasets  # <----- ⚠ this line is different to above script ⚠

# Following packages need to be independently installed. 
# ProgressMeter provides a progress bar while the user waits and Suppressor surpresses output to the terminal

@everywhere using ProgressMeter, Suppressor # <----- ⚠ this line is different to above script ⚠


# load data
tobs, yobs, σobs, _ = readdataset(source="Mrk6")

# We define an array of candidate delay vectors. Without loss of generalisation, the delay that corresponds to the first light curve is fixed to 0.
delays = [[0;d] for d in 0:0.2:100]

# We want to run cross-validation for all candidate delay vectors in parallel!

# Do "warmup" first for Julia
@showprogress pmap(D -> (@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, iterations=1000, numberofrestarts=3, delays = D, kernel = GPCC.matern32)), delays[1:2*nworkers()]) <----- ⚠ this line is different to above script, we use pmap instead map ⚠

# Do proper run 
out = @showprogress pmap(D -> (@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, iterations=2000, numberofrestarts=1, delays = D, kernel = GPCC.matern32)), delays) <----- ⚠ this line is different to above script, we use pmap instead map ⚠

# estimate posterior probability
getprobabilities(out)
```
