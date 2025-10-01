using Rasters, ArchGDAL

file_basename(year) = joinpath(
    @__DIR__, "../data/tif", "data_main_component_$year"
)

load_cycles(year) = DataFrame(Arrow.Table("$(file_basename(year)).cycles.arrow"))
load_summary(year) = DataFrame(Arrow.Table("$(file_basename(year)).summary.arrow"))
load_data(year) = Matrix{Float32}(Raster("$(file_basename(year)).tif").data)
