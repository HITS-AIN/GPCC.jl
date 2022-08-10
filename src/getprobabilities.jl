function getprobabilities(out)

    flatpriorvalues = ones(size(out))

    getprobabilities(out, flatpriorvalues)

end

function getprobabilities(out, priorpdfvalues)

    # add log prior to form log joint likelihood

    for i in eachindex(out)
        
        out[i] .+= log(priorpdfvalues[i])

    end

    Q = reduce(hcat, out)'

    nfolds = size(Q, 2)

    aux = vec(mean(reduce(hcat, [exp.(Q[:,i] .- logsumexp(Q[:,i])) for i in 1:nfolds]), dims=2))

    reshape(aux, size(out))
    
end
