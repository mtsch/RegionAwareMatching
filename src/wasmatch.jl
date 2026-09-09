include("data-loading.jl")

function birth_vertex(interval)
    # the matching code will create new intervals on the diagonal to match to. This replaces
    # them with missing
    if hasproperty(interval, :birth_simplex)
        return to_indices(interval.birth_simplex)
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
function main(; cutoff=1e-2)
    OUT_DIR = mkpath(joinpath(@__DIR__, "../data/match/"))
    years = 1990:2020
    death_cutoff = -0.05
    pers_cutoff = 0.01
    cutoff_func(x) = persistence(x) >= pers_cutoff && birth(x) <= death_cutoff
    for year in 1990:2019
        @info "$year"
        p1 = load_diagram(year)
        p2 = load_diagram(year+1)

        filter!(cutoff_func, p1.intervals)
        filter!(cutoff_func, p2.intervals)   
        for itv in p1.intervals
            filter!(cutoff_func, itv.children)
        end
        for itv in p2.intervals
            filter!(cutoff_func, itv.children)
        end

        output = matching_to_df(matching(Wasserstein(), p1, p2))
        Arrow.write("$OUT_DIR/match_pers$(cutoff)_$(year)_$(year+1).arrow", output)
    end
    @info "fin."
end

if !isinteractive()
    main()
end
