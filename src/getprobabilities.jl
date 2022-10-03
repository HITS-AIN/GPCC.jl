function getprobabilities(loglikel)

    flatpriorvalues = ones(size(loglikel))

    getprobabilities(loglikel, flatpriorvalues)

end

function getprobabilities(loglikel, logpriorpdfvalues)

    # add log prior to form log joint likelihood

    for i in eachindex(loglikel)
        
        loglikel[i] += logpriorpdfvalues[i]

    end

    posterior = exp.(loglikel .- logsumexp(loglikel))

    return posterior
    
end
