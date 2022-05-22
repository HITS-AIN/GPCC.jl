"""
out = performcv(tobs, yobs, σobs; delays = delays, iterations = 1, seedcv = 1, kernel = kernel, numberofrestarts = 1, numberoffolds = 5, plotting = false, rhomin = 0.1, rhomax = 20.0)

Performs cross-validation for Gaussian Process Cross Correlation (GPCC) model.

Data passed to the function are organised as arrays of arrays.
The outer array contains L number of inner arrays where L is the number of bands.
The l-th inner arrays hold the data pertaining to the l-th band.

See also [`getprobabilities`](@ref).

Input arguments
===============

- `tobs`: Array of arrays of observation times. There are L number of inner arrays. The l-th array holds the observation times of the l-th band.
- `yobs`: Array of arrays of fluxes. Same structure as `tobs`
- `σobs`: Array of error measurements. Same structure as `tobs`
- `kernel`: Specifies GP kernel function. Options are GPCC.OU / GPCC.rbf / GPCC.matern32
- `delays`: L-dimensional vector of delays.
- `iterations`: maximum number of iterations done when optimising marginal-likelihood of GP, i.e. optimising hyperparameters.
- `seedcv`: Random seed that controls the random sampling of initial solutions for GPCC and the random generation of cross-validation folds.
- `numberofrestarts`: Number of times to repeat optimisation in order to avoid suboptimal solutions due to poor initialisation.
- `initialrandom`: Before optimisation begins, a number of random solutions is sampled and the one with the highest likelihood becomes the starting point for the optimisation.
- `rhomin`: minimum value for lengthscale ρ of GP.
- `rhomax`: maximum value for lengthscale ρ of GP.

Output
===============

out: vector of test log-likelihoods for each fold

# Example
```julia-repl

julia> tobs, yobs, σobs = simulatedata(); # simulate data with default delays [0;2;6]
julia> out1 = performcv(tobs, yobs, σobs, iterations=1000, numberofrestarts=3, delays = [0;2;6], kernel = GPCC.matern32); # perform CV with true delays
julia> out2 = performcv(tobs, yobs, σobs, iterations=1000, numberofrestarts=3, delays = [0;2.1;5.9], kernel = GPCC.matern32); # julia> # perform CV with perturbed delays
julia> getprobabilities([out1, out2]) # estimate posterior probabilities, first entry corresponding to true delays should be higher
```
"""
function performcv(tobs, yobs, σobs; delays = delays, iterations = 1, seedcv = 1, kernel = kernel, numberofrestarts = 1, numberoffolds = 5, plotting = false, rhomin = 0.1, rhomax = 20.0)
    
    # let user know what is run
    str = @sprintf("\nRunning CV with %d number of folds\n\n", numberoffolds)
    colourprint(str, foreground = :cyan, bold = true)



    numberofbands = length(yobs)

    @assert(length(tobs) == length(yobs) == length(σobs))

    for i in 1:numberofbands
        @assert(length(tobs[i]) == length(yobs[i]) == length(σobs[i]))
    end


    # Specify folds for cross validation
    # Each band get its own partitioning

    bandfold = Array{Array{Array{Int64, 1},1}}(undef, numberofbands)

    for bandindex in 1:numberofbands
        rgkfold       = MersenneTwister(seedcv + bandindex)
        KF            = Kfold(length(tobs[bandindex]), numberoffolds)
        # problem with Kfold is that one cannot control the randomness of the partitioning
        # Hence the following workaround to what we would have normally used:
        # collect(Kfold(length(tobs[bandindex]), numberoffolds))
        KF.permseq   .= randperm(rgkfold, length(tobs[bandindex]))
        bandfold[bandindex] = collect(KF)
    end


    # Store here CV log-likelihoods calculated on left out fold
    fitness = zeros(numberoffolds)


    for foldindex in 1:numberoffolds


        # Specify training indices
        trainidx = [bandfold[bandindex][foldindex] for bandindex in 1:numberofbands]

        # Specify testing indices
        testidx = [sort(setdiff(1:length(tobs[bandindex]), trainidx[bandindex])) for bandindex in 1:numberofbands]

        # Sanity check: check splitting
        for i in 1:numberofbands
            @assert(all(length(trainidx[i]) + length(testidx[i]) == length(tobs[i])))
            @assert(Set(union(trainidx[i], testidx[i])) == Set(1:length(tobs[i])))
        end

        # split data into training and testing
        ttrain = [tobs[i][trainidx[i]] for i in 1:numberofbands]
        ytrain = [yobs[i][trainidx[i]] for i in 1:numberofbands]
        strain = [σobs[i][trainidx[i]] for i in 1:numberofbands]

        ttest  = [tobs[i][testidx[i]]  for i in 1:numberofbands]
        ytest  = [yobs[i][testidx[i]]  for i in 1:numberofbands]
        stest  = [σobs[i][testidx[i]]  for i in 1:numberofbands]

        str = @sprintf("\n--- fold %d train size is %d, test size is %d ---\n", foldindex, sum(length.(ttrain)), sum(length.(ttest)))
        colourprint(str, foreground = :light_blue, bold = true)

        # Run deconvolution
        predict =  @suppress gpcc(ttrain, ytrain, strain; kernel = kernel, delays = delays, numberofrestarts = numberofrestarts, iterations = iterations, seed = seedcv, rhomin = rhomin, rhomax = rhomax)[2]

        # evaluate on held out test data
        fitness[foldindex] = predict(ttest, ytest, stest)

        # produce plots for fold
        if plotting

            figure() ; cla()

            clr = ["r","g","b","m"]

            xtest = collect(minimum(map(minimum, tobs)):0.5:maximum(map(maximum, tobs)))

            local meanpred, Sigmapred = predict(xtest)

            for l in 1:numberofbands

                plot(ttrain[l], ytrain[l], "o"*clr[l])

                plot(ttest[l], ytest[l],   "x"*clr[l])


                local σpred = sqrt.(max.(Sigmapred[l], 1e-6))

                plot(xtest, meanpred[l], "--"*clr[l])

                fill_between(xtest, meanpred[l] .- σpred, meanpred[l] .+ σpred, color = clr[l], alpha=0.15)

            end

        end

        @printf("\t perf for fold %d is %f\n", foldindex, fitness[foldindex])

    end

    return fitness

end
