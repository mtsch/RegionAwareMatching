using Rasters, Ripserer, ArchGDAL, Plots
using ImageMorphology
using BenchmarkTools

function extract_main_component(year)
    raster_original = Raster("dhps_mapping_data/540_rasters/Resistance_median_$(year).tif"); 
    raster = copy(raster_original.data)
    
    # set all resistance values equal to one
    raster[raster .>= 0] .= 1
    # set all negative values to zero
    raster[raster.< 0] .= 0

    labelled_array = label_components(raster); # default uses C4 diamond shaped connectivity
    vec_array = component_lengths(labelled_array);

    new_vec_array = vec_array[begin+1:end];

    maxval, maxval_index = findmax(new_vec_array);

    # outputs the index of the second largest element 
    main_component_indices = findall(x->x==maxval_index, labelled_array);

    final_array = ones(size(raster)) .* -Inf
    final_array[main_component_indices] .= raster_original.data[main_component_indices]

    # overwrite original raster to obtain information needed to plot in correct orientation 
    raster_original.data .= final_array


    # save as tif file
    write("data_main_component_$(year).tif", raster_original)
end 

for year in 1990:2020
    extract_main_component(year)
    print(year)
end