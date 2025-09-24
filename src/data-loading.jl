# ====================== #
# Data loading shortcut. #
# ====================== #
using Rasters, ArchGDAL

file_basename(year) = joinpath(
    @__DIR__, "../data", "Resistance_median_main_component_$year"
)

"""
    load_data(year; raster=false)

Load data by year. If `raster` is set to `true`, return a `Raster` object, otherwise return
a `Matrix{Float32}`.

The data needs to be placed in the `data` directory for this to work.
"""
function load_data(year; raster=false)
    res = Raster("$(file_basename(year)).tif")
    if raster
        return res
    else
        return Matrix{Float32}(res.data)
    end
end
