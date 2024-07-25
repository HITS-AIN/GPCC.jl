function inducingpoints(x; dx = dx)

    xsorted = sort(x)

    # return xsorted[1]:dx:xsorted[end]

    close_enough(g) = minimum(abs.(g .- xsorted)) < dx
    
    grid = xsorted[1]:dx:xsorted[end]

    filter(close_enough, grid)

end
    