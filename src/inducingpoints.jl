function inducingpoints(x; dx = dx)

    z = Float64[]

    xsorted = sort(x)

    for xᵢ in xsorted

        if isempty(z)

            push!(z, xᵢ)

            continue

        end

        if abs(xᵢ-z[end]) > dx

            push!(z, xᵢ)

        end

    end

    return z

end