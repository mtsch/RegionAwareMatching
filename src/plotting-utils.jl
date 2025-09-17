include("shortest-rep.jl")

"""
    interval_str(interval)

Stringify the interval nicely for plotting.
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

"""
    plot_cycle!(ax, interval; threshold=death(interval), birth_simplex=true, kwargs...)

Plot an interval as a cycle at a given `threshold`. `kwargs` are passed to `lines!`.
"""
function plot_cycle!(
    ax, interval;
    linestyle=:solid,
    threshold=death(interval), birth_simplex=true, merge_tree=false, kwargs...
)
    cycle = minimum_area_cycle(interval; threshold)
    lines!(ax, cycle; linestyle, kwargs...)
    birth_vertex = Tuple(only(vertices(interval.birth_simplex)))
    if !isnothing(interval.parent)
        parent_vertex = Tuple(only(vertices(interval.parent_simplex)))
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

"""
    plot_heatmap!(fig, data; log=false, colormap=:viridis, kwargs...)

Prepare the base heatmap and plot it to `fig`.
"""
function plot_heatmap!(fig, data; log=false, colormap=:viridis, colorbar=true, kwargs...)
    heat = map(x -> !isfinite(x) ? missing : x, data)

    ax = Axis(fig[1, 1]; aspect=1, kwargs...)
    if log
        hm = image!(ax, log10.(heat); colormap)
        if colorbar
            Colorbar(fig[1, 2], hm; label=L"log10 value$$")
        end
    else
        hm = image!(ax, heat; colormap)
        if colorbar
            Colorbar(fig[1, 2], hm; label=L"value$$")
        end
    end

    ax.yreversed[] = true
    return ax
end

"""
    plot_merge_tree_leaf_segmentation(data, diagram; merge_tree=true, kwargs...)


Plot each interval in `diagram` and possibly a merge tree. Each cycle is plotted at a time
just before its first child merges into it.

Defined in `doi/10.1111/tgis.12816`.
"""
function plot_merge_tree_leaf_segmentation(diag; kwargs...)
    data = -diag.filtration.data
    plot_merge_tree_leaf_segmentation(data, diag; kwargs...)
    return fig
end
function plot_merge_tree_leaf_segmentation(data, diag; kwargs...)
    fig = Figure(size=(1920, 1080))
    plot_merge_tree_leaf_segmentation!(fig, data, diag; kwargs...)
    return fig
end
function plot_merge_tree_leaf_segmentation!(fig, diag;kwargs...)
    data = -diag.filtration.data
    plot_merge_tree_leaf_segmentation!(fig, data, diag; kwargs...)
end
function plot_merge_tree_leaf_segmentation!(fig, data, diag; merge_tree=true, legend=true, kwargs...)

    ax = plot_heatmap!(fig, data; kwargs...)

    for (i, interval) in enumerate(diag)
        color = Cycled(i ≥ 3 ? i + 1 : i)
        if isempty(interval.children)
            threshold = death(interval)
        else
            threshold = minimum(death, interval.children)
            @show threshold
        end
        label = L"$%$(interval_str(interval))$ at %$threshold"
        color = Cycled(i ≥ 3 ? i + 1 : i)

        plot_cycle!(ax, interval; threshold, label, color, merge_tree)

        if isfinite(death(interval)) && threshold ≠ death(interval)
            plot_cycle!(ax, interval; threshold=death(interval), label, color, merge_tree=false, linestyle=:dot)
        end
    end
    @show legend, merge_tree
    if legend && !isempty(diag)
        Legend(fig[:,3], ax; merge=true)
    end
    return ax
end

"""
    plot_all_at_threshold(data, diagram, threshold, kwargs...)

Plot each cycle at the same threshold.
"""
function plot_all_at_threshold(data, diagram, threshold; legend=true, kwargs...)
    fig = Figure(size=(1920, 1080))

    ax = plot_heatmap!(fig, data; kwargs...)

    for (i, interval) in enumerate(diagram)
        if birth(interval) > threshold
            continue
        end
        label = L"$%$(interval_str(interval))$"
        color = Cycled(i)
        plot_cycle!(ax, interval; threshold, color, label)
    end
    if legend && !isempty(diagram)
        Legend(fig[:,3], ax; merge=true)
    end
    return fig
end
