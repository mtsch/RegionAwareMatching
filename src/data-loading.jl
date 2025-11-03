# ====================== #
# Data loading shortcut. #
# ====================== #
using Rasters, ArchGDAL
using DataFrames
using Arrow

include("merge-trees.jl")

file_basename(year) = joinpath(
    @__DIR__, "../data/tif", "data_main_component_$year"
)

"""
    load_data(year; raster=false)

Load data by year. If `raster` is set to `true`, return a `Raster` object, otherwise return
a `Matrix{Float32}`.

The data needs to be placed in the `data` directory for this to work.
"""
function load_data(year; raster=false, lazy=true)
    res = Raster("$(file_basename(year)).tif"; lazy)
    if raster
        return res
    else
        return Matrix{Float32}(res.data)
    end
end

to_indices(::Nothing) = nothing
to_indices(c::Cube{0}) = Tuple(only(vertices(c)))
to_indices(c::Cube{1}) = Tuple.(vertices(c))
from_indices(_, ::Nothing) = nothing
from_indices(flt, indices::Tuple{Int,Int}) =
    simplex(flt, Val(0), (CartesianIndex(indices),))
from_indices(flt, indices::Tuple{Tuple,Tuple}) =
    simplex(flt, Val(1), CartesianIndex.(indices))

"""
    save_diagram(filename, diagram)

Save a diagram along with its merge tree information and representative to a compressed
arrow file.
"""
function save_diagram(filename, diagram)
    if isfinite(diagram.threshold)
        error("only infinite threshold supported")
    end

    index_map = Dict(zip(birth_simplex.(diagram), eachindex(diagram)))

    df = mapreduce(vcat, enumerate(diagram)) do (id, interval)
        representative = map(to_indices, interval.representative)
        parent = to_indices(interval.parent_simplex)
        birth = to_indices(interval.birth_simplex)
        death = to_indices(interval.death_simplex)

        DataFrame(; id, birth, death, parent, representative)
    end

    Arrow.write(filename, df; compress=:zstd)
end

"""
    load_diagram(diagram_file, data_file)
    load_diagram(year)

Load a diagram from arrow file created by `save_diagram` and data file, or by `year`.
"""
function load_diagram(year; lazy=true)
    diagram_file = joinpath(@__DIR__, "../data/diagrams/$year.arrow")
    data_file = "$(file_basename(year)).tif"

    return load_diagram(diagram_file, data_file; lazy)
end

function load_diagram(diagram_file, data_file; lazy=true)
    df = DataFrame(Arrow.Table(lazy ? diagram_file : read(diagram_file)))
    data = -Raster(data_file).data
    filtration = Cubical(data)

    intervals = PersistenceInterval[]
    for d in groupby(df, [:id])
        representative = map(d.representative) do vx
            simplex(filtration, Val(0), (vx,))
        end

        birth_simplex = from_indices(filtration, d.birth[1])
        death_simplex = from_indices(filtration, d.death[1])
        parent_simplex = from_indices(filtration, d.parent[1])

        interval = PersistenceInterval(
            birth(birth_simplex), isnothing(death_simplex) ? Inf : birth(death_simplex);
            birth_simplex,
            death_simplex,
            parent_simplex,
            children=PersistenceInterval[],
            representative,
        )
        push!(intervals, interval)
    end

    # Walk the merge tree and regenate its structure
    interval_map = Dict(zip(birth_simplex.(intervals), intervals))

    sort!(intervals; by=x -> (persistence(x), birth_simplex(x)), rev=true)

    # walk top-down adding parents to all intervals
    for i in eachindex(intervals)
        int = intervals[i]
        if !isnothing(int.parent_simplex)
            parent = interval_map[int.parent_simplex]
        else
            parent = nothing
        end
        int = PersistenceInterval(int.birth, int.death, (; int.meta..., parent))

        interval_map[int.birth_simplex] = int
        intervals[i] = int
    end

    # set children
    for int in intervals
        !isnothing(int.parent) && push!(int.parent.children, int)
    end
    sort!(intervals; by=x -> (persistence(x), birth_simplex(x)))

    return PersistenceDiagram(
        intervals; threshold=filtration.threshold, dim=0, field=Ripserer.Mod{2}, filtration
    )
end

# Test function - needs to be moved to testing suite
function check_identical(diagram1, diagram2)
    if diagram1.filtration.data ≠ diagram2.filtration.data
        error("filtration doesn't match")
    end
    diagram1 = sort(diagram1; by=x -> (persistence(x), birth_simplex(x)))
    diagram2 = sort(diagram2; by=x -> (persistence(x), birth_simplex(x)))
    if length(diagram1) ≠ length(diagram2)
        error("lengths don't match")
    end

    for i in eachindex(diagram1)
        left = diagram1[i]
        right = diagram2[i]
        for prop in [
            :birth, :death,
            :birth_simplex, :death_simplex,
            :parent_simplex, :parent,
            :representative,
        ]
            if getproperty(left, prop) ≠ getproperty(right, prop)
                error("property $prop at position $i does not match")
            end
        end
        if sort(left.children) ≠ sort(right.children)
            error("children at position $i don't match")
        end
    end
    return true
end

"""
    regenerate_diagrams()

Run persistent homology for all years and store the results to files in `data/diagrams`.
"""
function regenerate_diagrams()
    for year in 1990:2020
        @info "year $year"
        diagram = ripserer(
            Cubical(-load_data(year));
            merge_tree=true,
            dim_max=0,
            verbose=true,
            reps=true,
        )[1]
        save_diagram(joinpath(@__DIR__, "../data/diagrams/$year.arrow"), diagram)

        saved = load_diagram(year)
        @assert merge_tree_valid(saved)
        check_identical(diagram, saved)
    end
    @info "fin."
end
