include("data-loading.jl");
include("merge-trees.jl");
include("plotting-utils.jl");
include("tracking-utils.jl");

DATA_DIR = joinpath(@__DIR__, "../data/");

# original matching
# IMG_DIR = joinpath(@__DIR__, "../imgs/pers0.01");
# MATCH_DIR = joinpath(DATA_DIR, "match/match_pers0.01");
# MATCH_STR = "original matching";

# alternative matching
IMG_DIR = joinpath(@__DIR__, "../imgs/spatial_pers0.01");
MATCH_DIR = joinpath(DATA_DIR, "match/spatial_match_pers0.01");
MATCH_STR = "alternative matching";

LWD = 1.5

using ColorSchemes, Colors

# Palette 1: greedy selection
avoid_colors = sequential_palette(240, 500)[[100, 200, 300, 400]];
palette = distinguishable_colors(20, [RGB(0,0,0); RGB(1,1,1); avoid_colors], dropseed=true)[[1;3:end]]

# Palette 2: seaborn_colorblind
# palette = ColorSchemes.seaborn_colorblind
# palette = map(ColorSchemes.seaborn_colorblind) do color
#     hsl = HSL(color)
#     # Increase saturation, decrease lightness
#     new_saturation = min(hsl.s * 1.2, 1.0)
# 	new_lightness = hsl.l * 1.0
#     HSL(hsl.h, new_saturation, new_lightness)
# end
# palette = palette[2:9] # drop blues

# Palette 3: Set3_12
# palette = ColorSchemes.Set3_12.colors
# palette = map(ColorSchemes.Set3_12.colors) do color
#     hsl = HSL(color)
#     # Increase saturation, decrease lightness
#     new_saturation = min(hsl.s * 1.0, 1.0)
# 	new_lightness = hsl.l * 0.8
#     HSL(hsl.h, new_saturation, new_lightness)
# end
# palette = palette[[2:4;6:end]] # drop blues

set_theme!(palette=(color=palette,))

all_years = 1990:2020
pers_cutoff = 0.01

# Load persistence diagrams
@time diagrams = Dict(year => begin
    diagram = load_diagram(year; lazy=false)
    filter!(x -> persistence(x) >= pers_cutoff, diagram.intervals)
    for itv in diagram.intervals
        filter!(x -> persistence(x) >= pers_cutoff, itv.children)
    end
    diagram
end for year in all_years);

# Each dict maps birth vertex of interval to interval
vtx2itv_dicts = Dict(year => Dict(
	to_indices(itv.birth_simplex) => itv for itv in diagrams[year]
) for year in all_years);


# 1st tracking approach: track intervals (Vine)

vines = build_vines(all_years, vtx2itv_dicts, MATCH_DIR);
sort(countmap(length.(vines)))

# Plotting functions

function plot_vine(vine; display=true, save_dir=nothing, title_func=string)
    (save_dir !== nothing) && mkpath(save_dir)
    for (year, itv) in zip(vine.times, vine.itvs)
        raster = load_data(year; raster = true, lazy = false)
        fig = Figure()
        ax = plot_heatmap!(fig, raster; title=title_func(year))
        cycle, _ = minimum_area_cycle(itv)
        plot_cycle!(ax, raster, cycle; linewidth=LWD, color=palette[1])
        cycle, _ = minimum_area_cycle(itv.parent; threshold=itv.death)
        plot_cycle!(ax, raster, cycle; linewidth=LWD, linestyle=:dash, color=palette[1])

        display && display(fig)

        if save_dir !== nothing
            mkpath(save_dir)
            save(joinpath(save_dir, "$(year).png"), fig)
        end
        GC.gc()
    end
end

function plot_vines(vines; display=true, save_dir=nothing, title_func=string)
    yr_start = minimum((vine)->minimum(vine.times), vines)
    yr_end = maximum((vine)->maximum(vine.times), vines)
    (save_dir !== nothing) && mkpath(save_dir)
    for year in yr_start:yr_end
        raster = load_data(year; raster = true, lazy = false)
        fig = Figure()
        ax = plot_heatmap!(fig, raster; title=title_func(year))
        for (i, vine) in enumerate(vines)
            idx = findfirst(==(year), vine.times)
            idx === nothing && continue
            itv = vine.itvs[idx]
            cycle, _ = minimum_area_cycle(itv)
            plot_cycle!(ax, raster, cycle; linewidth=LWD, color=palette[i], alpha=0.8)
            cycle, _ = minimum_area_cycle(itv.parent; threshold=itv.death)
            plot_cycle!(ax, raster, cycle; linewidth=LWD, linestyle=:dash, color=palette[i], alpha=0.8)
        end
        display && display(fig)
        if save_dir !== nothing
            save(joinpath(save_dir, "$(year).png"), fig)
        end
        GC.gc()
    end
end

# Make series of plots corresponding to a vine

sort_by_totpers = sort(vines, by=(vine)->sum(persistence, vine.itvs), rev=true);
for idx in 2:3
    vine = sort_by_totpers[idx];
    plot_vine(
        vine; display=false, 
        save_dir=joinpath(IMG_DIR, "itv_totpers$idx"),
        title_func=(yr)->"Interval tracking ($yr, rank $idx by total persistence, $(MATCH_STR))"
    );
    GC.gc()
end
plot_vines(
    sort_by_totpers[2:5]; display=false, 
    save_dir=joinpath(IMG_DIR, "itv_totpers2-5"),
    title_func=(yr)->"Interval tracking ($yr, ranks 2-5 by total persistence, $(MATCH_STR))"
);


# 2nd tracking approach: track regions (RegionVine)

reg_vines = build_region_vines(all_years, vtx2itv_dicts, MATCH_DIR; info=true);
sort(countmap(length.(reg_vines)))

function plot_region_vines(vines; display=true, save_dir=nothing, title_func=string)
    yr_start = minimum((vine)->minimum(vine.times), vines)
    yr_end = maximum((vine)->maximum(vine.times), vines)
    (save_dir !== nothing) && mkpath(save_dir)
    for year in yr_start:yr_end
        raster = load_data(year; raster = true, lazy = false)
        fig = Figure()
        ax = plot_heatmap!(fig, raster; title=title_func(year))
        for (i, vine) in enumerate(vines)
            idx = findfirst(==(year), vine.times)
            idx === nothing && continue
            plot_cycle!(ax, raster, vine.cycles[idx]; linewidth=LWD, color=palette[i], alpha=0.8)
        end
        display && display(fig)
        if save_dir !== nothing
            save(joinpath(save_dir, "$(year).png"), fig)
        end
        GC.gc()
    end
end

sort_by_totpers = sort(reg_vines, by=(vine)->sum(persistence, vine.itvs), rev=true);
@assert all((x)->persistence(x)==Inf, sort_by_totpers[1].itvs)
plot_region_vines(
    sort_by_totpers[2:9]; display=false, 
    save_dir=joinpath(IMG_DIR, "reg_totpers2-9"),
    title_func=(yr)->"Region tracking ($yr, ranks 2-9 by total persistence, $(MATCH_STR))"
);

# Playground for investigating a specific vine
vine = sort_by_totpers[4];
yr1 = 1999; yr2 = yr1 + 1;
idx = findfirst(==(yr1), vine.times);
itv1, itv2 = vine.itvs[idx], vine.itvs[idx+1];
match_which(itv1, itv2)[1:3]

for year in (yr1, yr2)
    raster = load_data(year; raster = true, lazy = false)
    fig = Figure()
    ax = plot_heatmap!(fig, raster; title=string(year))
    idx = findfirst(==(year), vine.times)
    itv = vine.itvs[idx]
    cycle, _ = minimum_area_cycle(itv)
    plot_cycle!(ax, raster, cycle; linewidth=LWD, color=palette[1])
    cycle, _ = minimum_area_cycle(itv.parent; threshold=itv.death)
    plot_cycle!(ax, raster, cycle; linewidth=LWD, linestyle=:dash, color=palette[1])
    display(fig);
end