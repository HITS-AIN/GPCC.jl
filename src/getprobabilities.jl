function getprobabilities(out)

    Q = reduce(hcat, out)'

    nfolds = size(Q, 2)

    aux = vec(mean(reduce(hcat, [exp.(Q[:,i] .- logsumexp(Q[:,i])) for i in 1:nfolds]), dims=2))

    reshape(aux, size(out))
    
end
