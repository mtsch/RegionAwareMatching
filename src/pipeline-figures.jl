include("data-loading.jl");
include("merge-trees.jl");
include("plotting-utils.jl");
include("shortest-rep.jl");

all_years = 1990:2020
pers_cutoff = 0.01
death_cutoff = -0.05 # level sets will be clamped below death cutoff
cutoff_func(x) = persistence(x) >= pers_cutoff && birth(x) <= death_cutoff # birth should occur after death cutoff

DATA_DIR = joinpath(@__DIR__, "../data/");
IMG_DIR = joinpath(@__DIR__, "../imgs/pipeline_pers$(pers_cutoff)_death$(-death_cutoff)");
MATCH_DIR = joinpath(DATA_DIR, "match/region_match_pers$(pers_cutoff)_death$(-death_cutoff)");

mkpath(IMG_DIR)

year = 1990;
diagram = load_diagram(year; lazy=false);
filter!(cutoff_func, diagram.intervals); 
for itv in diagram.intervals
    filter!(cutoff_func, itv.children)
end

raster = load_data(year; raster = true, lazy = false);
segmentation = build_segmentation(diagram, raster; death_cutoff);

maximum(raster.data)
CMAX = 0.5

begin
    fig = Figure(size=(600, 500))
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
        title="Resistance prevalence in $year (>$(-death_cutoff) only)",
        titlesize=16
    )

    cmap_dict = Dict(itv => 0 for itv in diagram.intervals)
    hm = plot_segmentation!(
        ax, raster, segmentation, cmap_dict; 
        threshold=-death_cutoff, colorrange=(0, CMAX), domain_border=true
    )
    Colorbar(fig[1, 2], limits = (0, CMAX), colormap = GRAYMAP)
    display(fig)
    save(joinpath(IMG_DIR, "$(year)_masked_domain.png"), fig, px_per_unit=4)
end

itvs_by_area = sort(diagram.intervals, by=(itv)->length(segmentation[itv]), rev=true);
itvs_by_area[1:4]
itvs_by_area[2].death
itvs_by_area[3].death

thresholds = [0.2, 0.1, -itvs_by_area[3].death, -itvs_by_area[2].death, 0.05];
color_int_vec = [nothing, nothing, 3, 2, 1];

for (thres, color_int) in zip(thresholds, color_int_vec)
    thres_rounded = round(thres, digits=3)

    fig = Figure(size=(600, 500))
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
        title="Resistance prevalence in $year (>$(thres_rounded) only)",
        titlesize=16
    )

    cmap_dict = Dict(itv => 0 for itv in diagram.intervals)
    hm = plot_segmentation!(
        ax, raster, segmentation, cmap_dict; 
        threshold=thres, colorrange=(0, CMAX), domain_border=true
    )
    Colorbar(fig[1, 2], limits = (0, CMAX), colormap = GRAYMAP)
    display(fig)
    save(joinpath(IMG_DIR, "$(year)_bordered_$(thres_rounded).png"), fig, px_per_unit=4)
end

for (thres, color_int) in zip(thresholds, color_int_vec)
    thres_rounded = round(thres, digits=3)

    fig = Figure(size=(600, 500))
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
        title="Resistance prevalence in $year (>$(thres_rounded) only)",
        titlesize=16
    )

    cmap_dict = Dict(itv => 0 for itv in diagram.intervals)
    if !isnothing(color_int)
        cmap_dict[itvs_by_area[color_int]] = color_int
    end
    hm = plot_segmentation!(
        ax, raster, segmentation, cmap_dict; 
        threshold=thres, colorrange=(0, CMAX),
    )
    Colorbar(fig[1, 2], limits = (0, CMAX), colormap = GRAYMAP)
    display(fig)
    save(joinpath(IMG_DIR, "$(year)_levelset_$(thres_rounded).png"), fig, px_per_unit=4)
end


num_plot = 8
cmap_dict = Dict(itv => idx > 8 ? 0 : idx for (idx, itv) in enumerate(itvs_by_area));
begin
    fig = Figure(size=(500, 500));
    x_rad, y_rad = dims(raster)
    x_deg = rad2deg.(x_rad)
    y_deg = rad2deg.(y_rad)
    lyt = GridLayout(fig[1,1])
    ax = GeoAxis(
        lyt[1, 1];
        dest=EPSG(4326),
        limits=(extrema(x_deg), extrema(y_deg)),
        xticklabelsvisible=false, yticklabelsvisible=false,
        xgridvisible=false, ygridvisible=false,           
        title="Merge tree segmentation for $year",
    );
    plot_segmentation!(ax, raster, segmentation, cmap_dict)
    display(fig)
    save(joinpath(IMG_DIR, "$(year)_segmentation.png"), fig, px_per_unit=4)
end