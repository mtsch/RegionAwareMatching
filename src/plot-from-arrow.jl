using CairoMakie
using LaTeXStrings
using Arrow
using Rasters, ArchGDAL
using DataFramesMeta

file_basename(year) = joinpath(
    @__DIR__, "../data", "Resistance_median_main_component_$year"
)

load_cycles(year) = DataFrame(Arrow.Table("$(file_basename(year)).cycles.arrow"))
load_summary(year) = DataFrame(Arrow.Table("$(file_basename(year)).summary.arrow"))
load_data(year) = Matrix{Float32}(Raster("$(file_basename(year)).tif").data)

function interval_str((b, d))
    b = round(b, sigdigits=3)
    d = round(d, sigdigits=3)
    if isfinite(d)
        return L"[%$(b),%$(d))"
    else
        return L"[%$(b),\infty)"
    end
end

function plot_cycles(
    year;
    min_area=100, min_persistence=-Inf,
    dir=joinpath(@__DIR__, "../data"), padding=10, title=L"%$year",
    threshold=-0.02,
    colormap=:viridis,
    kwargs...
)
    heat = map(x -> !isfinite(x) ? missing : x, load_data(year))

    df = load_cycles(year)
    df = @rsubset df :area > min_area :threshold == threshold :persistence ≥ min_persistence
    sort!(df, [:area, :persistence]; rev=true)

    # Draw the heatmap
    fig = Figure(size=(1920, 1080))
    ax = Axis(fig[1, 1]; title, aspect=1, kwargs...)
    hm = heatmap!(ax, heat; colormap)
    Colorbar(fig[1, 2], hm; label=L"value$$")
    xlims!(ax, -padding, size(heat, 2) + padding)
    ylims!(ax, -padding, size(heat, 1) + padding)

    # Add cycles
    for (i, d) in enumerate(groupby(df, [:max_pos]))
        max_pos = d.max_pos[1]
        interval = d.interval[1]
        area = d.area[1]
        cycle = map(Point2, d.cycle)

        color = Cycled(i)
        linewidth = 2
        label = L"%$(interval_str(interval)), $A=%$area$"

        lines!(ax, cycle; label, linewidth, color)
        scatter!(ax, [max_pos]; label, color)
    end
    if !isempty(df)
        Legend(fig[:, 3], ax; merge=true)
    end

    xlims!(ax, (750,nothing))
    ylims!(ax, (400,nothing))
    ax.yreversed[] = true

    return fig
end

if false
    for year in 1990:2020
        for threshold in (-0.05, -0.5, -0.9)
            filename = joinpath(@__DIR__, "../plots", "$(year)_t$threshold.png")
            f = plot_cycles(year; threshold, min_area=1000)
            save(filename, f)
        end
    end
end
