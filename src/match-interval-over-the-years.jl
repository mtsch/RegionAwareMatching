## Scripts from the git repo
include("data-loading.jl")
include("merge-trees.jl")
include("plotting-utils.jl")

function select_first_child_of_largest(n=1; year=1990, cutoff=0.05)
    diagram = load_diagram(year)

    ## Sort by area, descending to make it easier to find the largest one
    sort!(diagram; by=area, rev=true)

    ## The diagram is a bit too big, so you should probably filter the tree. Below is the
    ## code that keeps the intervals that have persistence larger than 0.05.
    filtered = filter_merge_tree(x -> persistence(x) ≥ cutoff, diagram)

    largest = filtered[1]
    filtered = merge_subtree(filtered, largest)

    children = sort(largest.children; by=death, rev=true)
    return children[1:n]
end

function select_first_child_and_grandchild(; year=1990, cutoff=0.05)
    diagram = load_diagram(year)

    sort!(diagram; by=area, rev=true)
    filtered = filter_merge_tree(x -> persistence(x) ≥ cutoff, diagram)

    largest = filtered[1]
    filtered = merge_subtree(filtered, largest)

    child = argmax(death, largest.children)
    grandchild = argmax(death, child.children)

    return [child, grandchild]
end


function plot_a_cycle(year, intervals)
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
    for (i, interval) in enumerate(intervals)
        ismissing(interval) && continue
        cycle = map(minimum_area_cycle(interval)[1]) do (x, y)
            mean((x_deg[floor(Int,x)], x_deg[ceil(Int,x)])),
            mean((y_deg[floor(Int,y)], y_deg[ceil(Int,y)]))
        end

        plot!(ax, cycle, color=Cycled(i), markersize=2)
    end
    return fig
end

function track_interval(interval; file_postfix="_cut1.0e-5", start_at=1990)
    result = Vector{Union{Missing,PersistenceInterval}}(undef, length(start_at:2020))
    result .= missing
    interval_id = to_indices(interval.birth_simplex)

    result[1] = interval
    for (i, year) in enumerate(start_at:2020)
        year == start_at && continue

        matching = DataFrame(
            Arrow.Table("../data/match/match_$(year-1)_$(year)$(file_postfix).arrow"),
        )

        row = findfirst(matching.left) do id
            !ismissing(id) && id == interval_id
        end
        if isnothing(row)
            @warn "row not found in matching.left"
            return result
        end
        interval_id = matching[row, :right]

        if ismissing(interval_id)
            @info "Matched to the diagonal."
            return result
        end

        diagram = load_diagram(year)
        index_in_diagram = findfirst(diagram) do interval
            to_indices(interval.birth_simplex) == interval_id
        end
        interval = diagram[index_in_diagram]
        result[i] = interval
    end
    return result
end

function plot_matched_cycles(selected_interval::PersistenceInterval; kwargs...)
    return  plot_matched_cycles((selected_interval,); kwargs...)
end
function plot_matched_cycles(selected_intervals; name="wasmatch_", kwargs...)
    diagram = load_diagram(start_at)

    interval_sequences = map(x -> track_interval(x; kwargs...), selected_intervals)

    for (i, year) in enumerate(1990:2020)
        selected_intervals = map(x -> x[i], interval_sequences)
        fig = plot_a_cycle(year, selected_intervals)
        save("$(name)$(year).png", fig)
    end
end

if !isinteractive()
    # Old
    selected_interval = select_first_child_of_largest()
    plot_matched_cycles(selected_interval; name="one-child-")

    # Two ch's
    selected_intervals = select_first_child_of_largest(2; cutoff=0)
    plot_matched_cycles(selected_intervals; name="two-children-")

    #  Grandchildren
    selected_intervals = select_first_child_and_grandchild()
    plot_matched_cycles(selected_intervals; name="child-grandchild")
end
