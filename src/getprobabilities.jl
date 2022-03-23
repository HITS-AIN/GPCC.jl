function getprobabilities(out)

    Q = reduce(hcat, out)'

    nfolds = size(Q, 2)

    vec(mean(reduce(hcat, [exp.(Q[:,i] .- logsumexp(Q[:,i])) for i in 1:nfolds]), dims=2))

end
