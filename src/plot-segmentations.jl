include("region-wasmatch.jl");
include("plotting-utils.jl");

# Plot segmentations
for yr in 1990:2020
    diagram = load_diagram(yr; lazy=false);
    raster = load_data(yr; raster=true, lazy=false);

    pers_cutoff = 0.01
    cutoff_func(x) = persistence(x) >= pers_cutoff
    filter!(cutoff_func, diagram.intervals)  
    for itv in diagram.intervals
        filter!(cutoff_func, itv.children)
    end

    segmentation = build_segmentation(diagram, raster; death_cutoff=-0.05);
    itvs_by_area = sort(diagram.intervals, by=(itv)->length(segmentation[itv]), rev=true);

    num_plot = 8
    p = mod1.(3*yr .+ (1:num_plot), num_plot)
    cmap_dict = Dict(itv => idx > 8 ? 0 : p[idx] for (idx, itv) in enumerate(itvs_by_area))
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
            title="Merge tree leaf segmentation for $yr",
        );
        plot_segmentation!(ax, raster, segmentation, cmap_dict)
        # display(fig)
        path = joinpath(@__DIR__, "../imgs/segmentation_pers0.01_death0.05_area1-8")
        mkpath(path)
        save("$(path)/$(yr).png", fig, px_per_unit=4)
    end
end