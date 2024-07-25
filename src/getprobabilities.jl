"""
    pmf = getprobabilities(loglikel)

Converts log-likelihoods to probabilities.
"""
function getprobabilities(loglikel)

    flatpriorvalues = ones(size(loglikel))

    getprobabilities(loglikel, flatpriorvalues)

end

"""
    pmf = getprobabilities(loglikel, logpriorpdfvalues)

Converts log-likelihoods, accompanied by corresponding logpriorpdfvalues, to probabilities.
"""
function getprobabilities(loglikel, logpriorpdfvalues)

    # add log prior to form log joint likelihood

    joint_loglikel = loglikel .+ logpriorpdfvalues

    posterior = exp.(joint_loglikel .- logsumexp(joint_loglikel))

    return posterior
    
end
