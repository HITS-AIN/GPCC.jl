function Qvector(Narray, entries)

    L = length(Narray) # number of time series

    @assert(L == length(entries))

    local Q = zeros(sum(Narray))

    for l in 1:L

        Q[1+sum(Narray[1:l-1]):sum(Narray[1:l])] .= entries[l]

    end

    Q

end