function Qmatrix(Narray)

    L = length(Narray) # number of time series

    local Q = zeros(sum(Narray), L)

    for l in 1:L

        Q[1+sum(Narray[1:l-1]):sum(Narray[1:l]) ,l] .= 1.0

    end

    Q

end