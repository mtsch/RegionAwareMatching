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
        verbose=true,
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

function plot_a_cycle(year, interval)
    raster = Raster(joinpath(@__DIR__, "../data/Resistance_median_main_component_$(year).tif"))
    # replace -Inf with missing for nicer plot
    data = Matrix{Union{Missing,eltype(raster.data)}}(raster.data)
    data[.!isfinite.(data)] .= missing

    x_rad, y_rad = dims(raster)
    x_deg = rad2deg.(x_rad)
    y_deg = rad2deg.(y_rad)

    fig = Figure()
    ax = GeoAxis(
	fig[1, 1];
        dest=EPSG(4326),
	limits=(extrema(x_deg), extrema(y_deg)),
	xticklabelsvisible=false,
        xgridvisible=false,
	yticklabelsvisible=false,
        ygridvisible=false,
        title=L"%$year"
    )
    hm = heatmap!(
        ax,
	x_deg,
        y_deg,
        data,
	colorrange=(0, 1),
        colormap=sequential_palette(240, 500)[75:end]
    )
    lines!(GeoMakie.coastlines(), color="black")
    Colorbar(fig[1, 2], hm)

    # convert cycle to degrees
    cycle = map(minimum_area_cycle(interval)) do (x, y)
        mean((x_deg[floor(Int,x)], x_deg[ceil(Int,x)])),
        mean((y_deg[floor(Int,y)], y_deg[ceil(Int,y)]))
    end

    plot!(ax, cycle, color=:orangered3, markersize=2)
    return fig
end

function plot_matched_cycles(selected_interval)
    diagram = ripserer(
        Cubical(-load_data(1990); threshold=-0.05);
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

    fig = plot_a_cycle(1990, selected_interval)
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

        fig = plot_a_cycle(year, selected)
        save("wasmatch_$(year).png", fig)
    end
end

if !isinteractive()
    selected_interval = select_first_child_of_largest()
    plot_matched_cycles(selected_interval)
end
