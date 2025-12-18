include("data-loading.jl");
include("merge-trees.jl");
include("plotting-utils.jl");
include("shortest-rep.jl");

all_years = 1990:2020
pers_cutoff = 0.01
death_cutoff = -0.05 # level sets will be clamped below death cutoff
cutoff_func(x) = persistence(x) >= pers_cutoff && birth(x) <= death_cutoff # birth should occur after death cutoff

DATA_DIR = joinpath(@__DIR__, "../data/");
IMG_DIR = joinpath(@__DIR__, "../imgs/tracking_pers$(pers_cutoff)_death$(-death_cutoff)");
MATCH_DIR = joinpath(DATA_DIR, "match/region_match_pers$(pers_cutoff)_death$(-death_cutoff)");

# Functions

struct Vine{T<:Real}
	times::Vector{T}
	itvs::Vector{PersistenceInterval}
end
Base.length(vine::Vine) = length(vine.itvs)

function build_vines(years, vtx2itv_dicts, matching_prefix; info=false)
    # lookup uncompleted vines by most recent cycle, identified by birth vertex
    active = Dict{Tuple{Int64,Int64}, Vine}(); 
    completed = Vine[]; # completed vines

    for yr1 in years[1:end-1]
        yr2 = yr1 + 1
        to_extend = Dict{Tuple{Int64,Int64}, Vine}(); # vines that can continue to be to_extend
        matching_yr1 = DataFrame(
            Arrow.Table(read("$(matching_prefix)_$(yr1)_$(yr2).arrow")),
        )
        for (left, right) in eachrow(matching_yr1)
            (ismissing(left) || ismissing(right)) && continue
            itv1, itv2 = vtx2itv_dicts[yr1][left], vtx2itv_dicts[yr2][right]
            key1, key2 = only(vertices(itv1.birth_simplex)).I, only(vertices(itv2.birth_simplex)).I

            vine = pop!(active, key1, Vine([yr1], [vtx2itv_dicts[yr1][left]]))
            push!(vine.times, yr2)
            push!(vine.itvs, vtx2itv_dicts[yr2][right])
            to_extend[key2] = vine
        end
        append!(completed, collect(values(active)));
        if info
            @info "$yr1-$yr2" length(to_extend) length(completed)
        end
        active = to_extend
    end
    append!(completed, collect(values(active))); # complete all active vines
    return completed
end

function plot_vines(vines, segmentations; display_fig=true, save_dir=nothing, title_func=string, death_cutoff=Inf)
    yr_start = minimum((vine)->minimum(vine.times), vines)
    yr_end = maximum((vine)->maximum(vine.times), vines)
    (save_dir !== nothing) && mkpath(save_dir)
    for year in yr_start:yr_end
        raster = load_data(year; raster = true, lazy = false)
        fig = Figure()
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
            title=title_func(year),
        )

        diagram = diagrams[year]
        cmap_dict = Dict(itv => 0 for itv in diagram.intervals);
        for (i, vine) in enumerate(vines)
            idx = findfirst(==(year), vine.times)
            idx === nothing && continue
            cmap_dict[vine.itvs[idx]] = i            
        end

        segmentation = segmentations[year]
        plot_segmentation!(ax, raster, segmentation, cmap_dict)
        display_fig && display(fig)
        if save_dir !== nothing
            save(joinpath(save_dir, "$(year).png"), fig)
        end
        GC.gc()
    end
end

function plot_vine_persistence!(ax, vine, colormap)
    years = vine.times
    births = birth.(vine.itvs)
    deaths = min.(death.(vine.itvs), 0.)
    nc = length(colormap)
    band!(
        ax, years, -deaths, -births, 
        strokecolor=colormap[end], strokewidth=2, 
        color=(colormap[nc÷3], 0.3)
    )
end

function get_area(itv, year, segmentations)
    length(segmentations[year][itv])
end

function plot_vine_area!(ax, vine, segmentations; kwargs...)
    areas = [get_area(itv, year, segmentations) for (year, itv) in zip(vine.times, vine.itvs)]
    lines!(
        ax, vine.times, areas; kwargs...
    )
end

# Load persistence diagrams
@time diagrams = Dict(year => begin
    diagram = load_diagram(year; lazy=false)
    filter!(cutoff_func, diagram.intervals)   
    for itv in diagram.intervals
        filter!(cutoff_func, itv.children)
    end
    diagram
end for year in all_years);

# Compute all segmentations
segmentations = Dict(
    year => begin 
        diagram = diagrams[year]
        raster = load_data(year; raster = true, lazy = false)
        build_segmentation(diagram, raster; death_cutoff);
    end for year in all_years
);

# Each dict maps birth vertex of interval to interval
vtx2itv_dicts = Dict(year => Dict(
	to_indices(itv.birth_simplex) => itv for itv in diagrams[year]
) for year in all_years);

# Build and rank vines
vines = build_vines(all_years, vtx2itv_dicts, MATCH_DIR);
sort(countmap(length.(vines))) # summarise vine length
sort_by_totarea = sort(vines, by=(vine)->sum(get_area.(vine.itvs, vine.times, Ref(segmentations))), rev=true);

begin
    fig = Figure(size=(1200, 400))
    ax = Axis(
        fig[1,1], limits=(extrema(all_years), (0., 1.)), 
        xlabel="Year", ylabel="Drug resistance prevalence"
    )
    for k in 1:8
        vine = sort_by_totarea[k]
        plot_vine_persistence!(ax, vine, COLORMAPS[k])
    end
    display(fig)
    save(joinpath(IMG_DIR, "totarea1-8_persistences.png"), fig)
end

begin
    fig = Figure(size=(1200, 400))
    ax = Axis(fig[1,1], xlabel="Year", ylabel="Area (number of pixels)")
    for k in 1:8
        vine = sort_by_totarea[k]
        plot_vine_area!(ax, vine, segmentations; color=COLORMAPS[k][end])
    end
    display(fig)
    save(joinpath(IMG_DIR, "totarea1-8_areas.png"), fig)
end

plot_vines(
    sort_by_totarea[1:8], segmentations; death_cutoff, display_fig=true, 
    save_dir=joinpath(IMG_DIR, "totarea1-8"),
    title_func=(yr)->"Region-aware tracking ($yr, ranks 1-8 by total area)"
);
