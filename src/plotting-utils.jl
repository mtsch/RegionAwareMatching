using CairoMakie
using GeoMakie
using LaTeXStrings

include("shortest-rep.jl")

"""
    interval_str(interval)

Stringify the interval to LaTeX for plotting.
"""
function interval_str((b, d))
    b = round(b, sigdigits=3)
    d = round(d, sigdigits=3)
    if isfinite(d)
        return L"[%$(b),%$(d))"
    else
        return L"[%$(b),\infty)"
    end
end

function plot_heatmap!(
    ax::GeoAxis, raster::Raster;
    threshold=Inf,
    colormap=sequential_palette(240, 500)[75:end],
    bg_colormap=sequential_palette(240, 1000; s=0)[250:500],
    coastlines=true,
    kwargs...
)
    x_rad, y_rad = dims(raster)
    x_deg = rad2deg.(x_rad)
    y_deg = rad2deg.(y_rad)

    x = extrema(x_deg)
    y = reverse(extrema(y_deg))

    # replace -Inf with missing for nicer plot
    data = Matrix{Union{Missing,eltype(raster.data)}}(raster.data)
    data[.!isfinite.(data)] .= missing

    if threshold ≠ Inf

        # apply thresholding
        mask = isfinite.(data)
        data_masked = copy(data)
        data_mask = copy(data)
        data_masked[ismissing.(data) .|| threshold .> data] .= missing
        data_mask[ismissing.(data) .|| threshold .≤ data] .= missing

        image!(ax, x, y, data_mask; colormap=bg_colormap)
    else
        data_masked = data
    end


    hm = image!(ax, x, y, data_masked; colorrange=(0, 1), colormap)
    if coastlines
        lines!(GeoMakie.coastlines(), color="black")
    end

    return hm
end
function plot_heatmap!(fig, raster::Raster; title="", colorbar=true, kwargs...)
    x_rad, y_rad = dims(raster)
    x_deg = rad2deg.(x_rad)
    y_deg = rad2deg.(y_rad)

    lyt = GridLayout(fig[1,1])

    ax = GeoAxis(
        lyt[1, 1];
        dest=EPSG(4326),
        limits=(extrema(x_deg), extrema(y_deg)),
        xticklabelsvisible=false,
        xgridvisible=false,
        yticklabelsvisible=false,
        ygridvisible=false,
        title,
    )
    hm = plot_heatmap!(ax, raster; kwargs...)
    if colorbar
        Colorbar(lyt[1, 2], hm)
    end

    return ax
end
function plot_heatmap(raster::Raster; kwargs...)
    fig = Figure()
    plot_heatmap!(fig, raster; kwargs...)
    return fig
end

"""
    plot_cycle!(ax, data, interval; threshold=death(interval), birth_simplex=true, kwargs...)

Plot an interval as a cycle at a given `threshold`. `kwargs` are passed to `lines!`.
"""
function plot_cycle!(
    ax, data::Raster, interval;
    linestyle=:solid,
    linewidth=2,
    threshold=death(interval),
    birth_simplex=true,
    merge_tree=false,
    kwargs...
)
    x_rad, y_rad = dims(raster)
    x_deg = rad2deg.(x_rad)
    y_deg = rad2deg.(y_rad)

    cycle = minimum_area_cycle(interval; threshold=-threshold)

    remapped_cycle = map(cycle) do (x, y)
        (
            mean((x_deg[floor(Int, x)], x_deg[ceil(Int, x)])),
            mean((y_deg[floor(Int, y)], y_deg[ceil(Int, y)])),
        )
    end

    lines!(ax, remapped_cycle; linestyle, linewidth, kwargs...)
    birth_x, birth_y = Tuple(only(vertices(interval.birth_simplex)))
    birth_vertex = (x_deg[birth_x], y_deg[birth_y])

    if !isnothing(interval.parent)
        parent_x, parent_y = Tuple(only(vertices(interval.parent_simplex)))
        parent_vertex = (x_deg[parent_x], y_deg[parent_y])
    else
        parent_vertex = nothing
    end
    if birth_simplex
        scatter!(ax, [birth_vertex]; kwargs...)
    end
    if merge_tree && !isnothing(parent_vertex)
        lines!(ax, [birth_vertex, parent_vertex]; kwargs...)
    end
    return ax
end
function plot_cycle(data::Raster, interval; kwargs...)
    fig = Figure()
    ax = plot_heatmap!(fig, data)
    plot_cycle!(ax, data, interval; kwargs...)
    return fig
end

"""
    plot_merge_tree_leaf_segmentation(data, diagram; merge_tree=true, kwargs...)


Plot each interval in `diagram` and possibly a merge tree. Each cycle is plotted at a time
just before its first child merges into it.

Defined in `doi/10.1111/tgis.12816`.
"""
function plot_merge_tree_leaf_segmentation(data::Raster, diag; kwargs...)
    fig = Figure(size=(1920, 1080))
    plot_merge_tree_leaf_segmentation!(fig, data, diag; kwargs...)
    return fig
end
function plot_merge_tree_leaf_segmentation!(
    fig, data::Raster, diag; merge_tree=true, legend=true, kwargs...
)

    ax = plot_heatmap!(fig, data; kwargs...)

    for (i, interval) in enumerate(diag)
        color = Cycled(i ≥ 3 ? i + 1 : i)
        if isempty(interval.children)
            threshold = death(interval)
        else
            threshold = minimum(death, interval.children)
        end
        label = L"$%$(interval_str(interval))$ at %$threshold"
        color = Cycled(i ≥ 3 ? i + 1 : i)

        plot_cycle!(ax, data, interval; threshold, label, color, merge_tree)

        if isfinite(death(interval)) && threshold ≠ death(interval)
            plot_cycle!(ax, data, interval; threshold=death(interval), label, color, merge_tree=false, linestyle=:dot)
        end
    end
    if legend && !isempty(diag)
        Legend(fig[:,3], ax; merge=true)
    end
    return ax
end

"""
    plot_all_at_threshold(data, diagram, threshold, kwargs...)

Plot each cycle at the same threshold.
"""
function plot_all_at_threshold(
    data, diagram, threshold; legend=true, birth_simplex=false,
    xlims=nothing, ylims=nothing,
    colorbar=true,
    kwargs...
)
    fig = Figure()

    ax = plot_heatmap!(fig, data; threshold, colorbar, kwargs...)

    for (i, interval) in enumerate(diagram)
        if -death(interval) ≤ threshold ≤ -birth(interval)
            label = L"$%$(interval_str(interval))$"
            color = Cycled(i)
            plot_cycle!(ax, data, interval; threshold, color, label, birth_simplex)
        end
    end
    if legend && !isempty(diagram)
        Legend(fig[:, 2], ax; merge=true)
    end

    !isnothing(xlims) && xlims!(ax, xlims)
    !isnothing(ylims) && ylims!(ax, ylims)
    return fig
end

# ylims=(-20, 20)
# xlims=(10, 40)
