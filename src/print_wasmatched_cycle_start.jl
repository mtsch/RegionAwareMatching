## This code finds the cycle in year 1990 that we want to follow until 2020
### !!! When switching between julia and python be aware of index shift !!!
using Ripserer
using CSV

## Scripts from the git repo
include("data-loading.jl")
include("merge-trees.jl")
include("plotting-utils.jl")

function select_first_child_of_largest()
    data = load_data(1990)
    diagram = ripserer(
        Cubical(-data; threshold=-0.05);
        merge_tree=true,
        dim_max=0,
        reps=true,
    )[1]

    ## Sort by area, descending to make it easier to find the largest one
    sort!(diagram; by=area, rev=true)

    ## The diagram is a bit too big, so you should probably filter the tree. Below is the
    ## code that keeps the intervals that have persistence larger than 0.05.
    filtered = filter_merge_tree(x -> persistence(x) ≥ 0.05, diagram)

    largest = filtered[1]
    filtered = merge_subtree(filtered, largest)

    selected = argmin(death, largest.children)
    return selected
end

function plot_a_cycle(data, interval; kwargs...)
    fig = Figure(size=(800, 600))

    ax = plot_heatmap!(fig, data; kwargs...)
    plot!(ax, minimum_area_cycle(selected), color=:orangered3, markersize=2)

    return fig
end

function plot_matched_cycles(selected_interval)
    diagram = ripserer(
        Cubical(-data; threshold=-0.05);
        # merge_tree=true,
        dim_max=0,
        reps=true,
        verbose=true,
    )[1]

    selected_index = findfirst(diagram) do int
        int.birth_simplex == selected_interval.birth_simplex
    end
    @info "Year 1990 selected index: $selected_index"

    if isnothing(selected_index)
        error("interval $selected_interval not found in persistence diagram!")
    end

    fig = plot_a_cycle(data, selected_interval; title=L"1990")
    save("wasmatch_1990.png", fig)

    for year in 1991:2020
        matching = CSV.read("../data/matching_$(year-1)_$(year).csv", DataFrame)

        matching_dict = Dict{Int,Int}(zip(matching[:,1] .+ 1, matching[:,2] .+ 1))

        selected_index = matching_dict[selected_index]
        @info "Year $year selected index: $selected_index"

        diagram = first(
            ripserer(
                Cubical(-load_data(year); threshold=-0.05);
                merge_tree=true, dim_max=0, reps=true, verbose=true,
            )
        )
        selected = diagram[selected_index]

        fig = plot_a_cycle(data, selected; title=L"%$year")
        save("wasmatch_$(year).png", fig)
    end
end

if !isinteractive()
    selected_interval = select_first_child_of_largest()
    plot_matched_cycles(selected_interval)
end
