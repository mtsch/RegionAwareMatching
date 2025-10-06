include("data-loading.jl")

function birth_vertex(interval)
    # the matching code will create new intervals on the diagonal to match to. This replaces
    # them with missing
    if hasproperty(interval, :birth_simplex)
        return Tuple(only(vertices(interval.birth_simplex)))
    else
        return missing
    end
end

function matching_to_df(m)
    result = DataFrame()
    for (left, right) in matching(m)
        push!(result, (left=birth_vertex(left), right=birth_vertex(right)); promote=true)
    end
    return result
end

## produces matchings for all pairs of consecutive years
function main()
    for year in 1990:2020
        p1 = load_diagram(year)
        p2 = load_diagram(year+1)
        output = matching_to_df(matching(Wasserstein(), p1, p2))
        Arrow.write("match_$year_$(year+1).arrow", output)
    end
end
