## This code finds the cycle in year 1990 that we want to follow until 2020
### !!! When switching between julia and python be aware of index shift !!!

using Ripserer
using CSV

## Scripts from the git repo
include("data-loading.jl")
include("merge-trees.jl")
include("plotting-utils.jl")

data = load_data(1990)
diagram = first(ripserer(Cubical(-data; threshold=-0.05); merge_tree=true, dim_max=0, reps=true))

## Sort by area, descending to make it easier to find the largest one
sort!(diagram; by=area, rev=true)

## The diagram is a bit too big, so you should probably filter the tree. Below is the code
## that keeps the intervals that have persistence larger than 0.05.
filtered = filter_merge_tree(x -> persistence(x) ≥ 0.05, diagram)

## If you want you can also drop all intervals that are not in a subtree with the largest one
largest = filtered[1]
filtered = merge_subtree(filtered, largest)

## Now you can look at the children to find the one that merged into it first or last
## If we look at the min2, we get the child that is the first one to merge
println("Min Child:", minimum(death, largest.children))
println("Max Child:", maximum(death, largest.children))
println("Min Child's Child:", minimum(death, largest.children[1].children))
println("Max Child's Child:", maximum(death, largest.children[1].children))

## We want to find the cycle with the corresponding death times


selected = argmin(death, largest.children)

diagram = first(ripserer(Cubical(-data; threshold=-0.05); merge_tree=true, dim_max=0, reps=true, verbose=true))

selected_index = findfirst(diagram) do int
    int.birth_simplex == selected.birth_simplex
end

if isnothing(selected_index)
    error("you messed up again!")
end

fig = Figure(size=(800, 600))

plot_heatmap!(fig, data, title="1990")
plot!(minimum_area_cycle(selected), color=:orangered3, markersize=2)

save("wasmatch_1990.png", fig)

@show selected_index
for year in 1991:2020
    matching = CSV.read("../data/matching_$(year-1)_$(year).csv", DataFrame)

    matching_dict = Dict{Int,Int}(zip(matching[:,1] .+ 1, matching[:,2] .+ 1))

    selected_index = matching_dict[selected_index]

    diagram = first(
        ripserer(
            Cubical(-load_data(year); threshold=-0.05);
            merge_tree=true, dim_max=0, reps=true, verbose=true,
        )
    )
    selected = diagram[selected_index]

    fig = Figure(size=(800, 600))

    ax = plot_heatmap!(fig, data, title="$year")
    plot!(minimum_area_cycle(selected), color=:orangered3, markersize=2)

    save("wasmatch_$(year).png", fig)

    @show selected_index
end
