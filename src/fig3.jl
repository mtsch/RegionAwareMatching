using Ripserer

include("data-loading.jl")
include("merge-trees.jl")
include("plotting-utils.jl")



function pick_top_3_cycles(data; min_persistence=0.05)
    diag = first(ripserer(
        Cubical(-data; threshold=-min_persistence);
        dim_max=0, merge_tree=true, verbose=true, reps=true
    ))
    largest = argmax(area, diag)
    filtered = merge_subtree(diag, largest)
    filter_merge_tree!(x -> persistence(x) ≥ min_persistence, filtered)
    return sort!(diag; by=area, rev=true)
    #@assert merge_tree_valid(filtered)
    #return sort!(filtered; by=area, rev=true)
end

function plot_at_time!(ax, data, diagram, threshold)
    for (i, int) in enumerate(diagram)
        if birth(int) ≤ threshold ≤ death(int)
            plot_cycle!(ax, data, int; threshold, color=Cycled(i ≥ 3 ? i + 1 : i), birth_simplex=false)
        elseif death(int) ≤ threshold
            plot_cycle!(ax, data, int; threshold, color=Cycled(i ≥ 3 ? i + 1 : i), linestyle=:dot, birth_simplex=false)
        end
    end
end

function big_plot(year; kwargs...)
    raster = load_data(year; raster=true)
    data = load_data(year)
    diagram = pick_top_3_cycles(data)
    return big_plot(raster, diagram; kwargs...)
end

function big_plot(data, diagram; xlims=(900, nothing), ylims=(400, nothing))
    sort!(diagram, by=area, rev=true)
    top3 = deepcopy(merge_subtree(diagram, diagram[1]))
    sort!(top3, by=area, rev=true)
    selected_area = area(top3[3])
    filter_merge_tree!(x -> area(x) ≥ selected_area, top3)
    sort!(top3, by=area, rev=true)
    display(top3)
    top3 = top3[1:3]

    fig = Figure(size=(800, 600))

    fig[1,1] = subgl_top = GridLayout()

    # Plot snapshots
    kws = (;xticklabelsvisible=false, yticklabelsvisible=false)
    top_axes = Vector{GeoAxis}(undef, 4)
    top_axes[1] = plot_heatmap!(subgl_top[2,1], data; colorbar=false, kws...)
    top_axes[2] = plot_heatmap!(subgl_top[2,2], data; colorbar=false, kws...)
    top_axes[3] = plot_heatmap!(subgl_top[2,3], data; colorbar=false, kws...)
    top_axes[4] = plot_heatmap!(subgl_top[2,4], data; colorbar=false, kws...)

    times = reverse!(sort([
        0.05, #=0.5,=# 0.9, -death(top3[2]), -death(top3[3])
    ]))

    for (i, (ax, time)) in enumerate(zip(top_axes, times))
        Label(subgl_top[1,i], L"t=%$(round(time; digits=3))")
        plot_at_time!(ax, data, top3, -time)
        xlims!(ax, xlims...)
        ylims!(ax, ylims...)
        ax.yreversed[] = true
        colsize!(subgl_top, i, Relative(1/6))
    end

    fig[2:3,1] = subgl_bot = GridLayout()

    # MTLS
    ax = plot_merge_tree_leaf_segmentation!(subgl_bot[1,2:4], top3; legend=false, merge_tree=false)
    Box(subgl_bot[1,1]; color=RGBAf(0,0,0,0), strokevisible=false)
    Box(subgl_bot[1,5]; color=RGBAf(0,0,0,0), strokevisible=false)
    xlims!(ax, xlims...)
    ylims!(ax, ylims...)
    ax.yreversed[] = true

    # diagram
    # plot_diagram!(fig[3,1:5], diagram)

    return fig
end
