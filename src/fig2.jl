include("data-loading.jl")
include("plotting-utils.jl")

function diagram_plot(diagram, selected)
    diagram = filter(diagram) do int
        int.birth_simplex ∉ birth_simplex.(selected)
    end

    fig = Figure(size=(370, 300))

    hi = -minimum(b for (b, _) in diagram)
    lo = -maximum(d for (_, d) in diagram if isfinite(d))
    inf = round(4lo) / 4

    xticks = (0:0.2:1.2, [i == inf ? L"-∞" : "$i" for i in 0:0.2:1.2])

    ax = Axis(fig[1, 1]; xlabel=L"death$$", ylabel=L"birth$$", aspect=1, xticks)

    lines!(ax, [inf,hi], [inf,hi]; linestyle=:dash, color=:black)
    vlines!(ax, [inf]; linestyle=:dot, color=:gray)

    diag_points = [Point2f(isfinite(d) ? -d : inf, -b) for (b, d) in diagram]
    sel_points = [Point2f(isfinite(d) ? -d : inf, -b) for (b, d) in selected]

    scatter!(ax, diag_points; color=:gray, markersize=2, marker=:diamond)
    for (i, pt) in enumerate(sel_points)
        color = Cycled(i + 1)
        scatter!(ax, [pt]; color)
    end

    return fig
end

function prepare_fig2(
    diagram=nothing; year=2002, thresholds=(0.15, 0.5, 0.7, 0.9), min_area=nothing
)

    raster = load_data(year; raster=true)

    if isnothing(diagram)
        diagram = ripserer(
            Cubical(-raster.data; threshold=-0.05);
            dim_max=0, reps=true, verbose=true, merge_tree=true,
        )[1]
    end

    # Get main connected component
    selected = merge_subtree(diagram, argmax(area, diagram))

    # Only keep intervals that will appear on the plot
    selected = filter(diagram) do int
        any(thresholds) do t
            -death(int) ≤ t ≤ -birth(int)
        end
    end
    sort!(selected, by=area; rev=true)
    selected = filter(diagram) do int
        if isnothing(min_area)
            # only keep top 3
            area(int) ≥ area(selected[3])
        else
            area(int) > min_area
        end
    end
    display(collect(zip(selected, area.(selected))))

    for threshold in thresholds
        @show threshold
        f = plot_all_at_threshold(
            raster, selected, threshold;
            colorbar=false, legend=false, #ylims=(-20, 20), xlims=(10, 40),
        )
        save("$(year)_$(threshold).pdf", f)
    end
    f = diagram_plot(diagram, selected)
    save("$(year)_diagram.pdf", f)
end
