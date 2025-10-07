## Scripts from the git repo
include("data-loading.jl")
include("merge-trees.jl")
include("plotting-utils.jl")

function select_first_child_of_largest(year=1990)
    diagram = load_diagram(year)

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
    raster = load_data(year; raster=true)
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

function plot_matched_cycles(selected_interval; file_postfix="_cut1.0e-5", start_at=1990)
    diagram = load_diagram(start_at)

    fig = plot_a_cycle(start_at, selected_interval)
    save("wasmatch_$(start_at).png", fig)

    selected_interval_id = to_indices(selected_interval.birth_simplex)

    for year in start_at+1:2020
        matching = DataFrame(
            Arrow.Table("../data/match/match_$(year-1)_$(year)$(file_postfix).arrow"),
        )

        row = findfirst(matching.right) do id
            !ismissing(id) && id == selected_interval_id
        end
        selected_interval_id = matching[row, :right]

        if ismissing(selected_interval_id)
            @info "Matched to the diagonal."
            return
        end


        diagram = load_diagram(year)
        index_in_diagram = findfirst(diagram) do interval
            to_indices(interval.birth_simplex) == selected_interval_id
        end
        selected = diagram[index_in_diagram]

        @info "Year $year selected id: $selected_interval_id ($(selected))"

        fig = plot_a_cycle(year, selected)
        save("wasmatch_$(year).png", fig)
    end
end

if !isinteractive()
    selected_interval = select_first_child_of_largest()
    plot_matched_cycles(selected_interval)
end
