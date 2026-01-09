# ============================================================= #
# Code for extracting minimum area cycles from representatives. #
# ============================================================= #
using Ripserer, PersistenceDiagrams
using DelaunayTriangulation
using StaticArrays
using DataStructures
using DataFramesMeta
using Arrow
using Rasters, ArchGDAL, JSON

# 4-connectivity offsets
OFFSETS = (
    CartesianIndex(1, 0),
    CartesianIndex(-1, 0),
    CartesianIndex(0, 1),
    CartesianIndex(0, -1)
)

# Merge tree leaf segmentation for region-aware tracking (https://arxiv.org/abs/2510.16486)

SegDict = Dict{PersistenceInterval,Set{CartesianIndex{2}}};

function build_segmentation(diagram::PersistenceDiagram, raster::Raster; death_cutoff=Inf)::SegDict
    itvs_by_birth = sort(diagram.intervals, by=birth)
    idx_by_itv = Dict(reverse.(enumerate(itvs_by_birth)))
    segment_arr = zeros(Int, size(raster.data));
    pqueue = PriorityQueue();
    for (idx, itv) in enumerate(itvs_by_birth)
        vertex = only(vertices(itv.birth_simplex))
        enqueue!(pqueue, (vertex, idx), (-raster[vertex], idx))
    end
    while !isempty(pqueue)
        curr_vertex, curr_idx = dequeue!(pqueue)
        segment_arr[curr_vertex] == 0 || continue # skip if vertex is already assigned a region
        segment_arr[curr_vertex] = curr_idx
        for offset in OFFSETS
            neighbour = curr_vertex + offset
            val = -raster[neighbour]
            if val < death_cutoff && segment_arr[neighbour] == 0
                idx = val < death(itvs_by_birth[curr_idx]) ? curr_idx : idx_by_itv[itvs_by_birth[curr_idx].parent]
                pqueue[(neighbour, idx)] = (val, idx)
            end
        end
    end

    segmentation = Dict(itv => Set{CartesianIndex}() for itv in itvs_by_birth);
    for I in CartesianIndices(segment_arr)
        if segment_arr[I] > 0
            push!(segmentation[itvs_by_birth[segment_arr[I]]], I)
        end
    end

    return segmentation
end

"""
    threshold_representative(filtration, interval, thresh)

Prune the representative of the interval to only include vertices that are present at a
given threshold. Returns a `PersistenceInterval` with the new representative.
"""
function threshold_representative(filtration::Ripserer.AbstractFiltration, int, thresh)
    rep_vertices = threshold_representative(int.representative, int.birth_simplex, thresh)
    representative = map(v -> simplex(filtration, Val(0), (v,)), rep_vertices)
    return PersistenceInterval(int.birth, int.death, (; int.meta..., representative))
end
function threshold_representative(representative::AbstractVector, birth_simplex, thresh)
    rep_vertices = Set(only(vertices(s)) for s in representative if birth(s) < thresh)
    visited = empty(rep_vertices)
    stack = [only(vertices(birth_simplex))]
    representative = eltype(rep_vertices)[]

    while !isempty(stack)
        curr = pop!(stack)
        push!(representative, curr)
        for dir in ((-1,0), (1,0), (0,-1), (0,1))
            neighbour = curr + CartesianIndex(dir)
            if neighbour ∈ rep_vertices && neighbour ∉ visited
                push!(visited, neighbour)
                push!(stack, neighbour)
            end
        end
    end
    return representative
end

function representative_mask(interval::PersistenceInterval, threshold=Inf)
    return representative_mask(interval.representative, interval.birth_simplex, threshold)
end
function representative_mask(representative, birth_simplex, threshold=Inf)
    rep = threshold_representative(representative, birth_simplex, threshold)

    height = maximum(v -> v[1], rep)
    width = maximum(v -> v[2], rep)
    mask = zeros(Bool, (height, width))
    for c in rep
        mask[c] = true
    end
    return mask
end

struct Neighbourhood <: AbstractMatrix{Bool}
    mask::SMatrix{3,3,Bool}
end
function Neighbourhood(mask, index)
    hood = @SMatrix[get(mask, index + CartesianIndex(i, j), false) for i in -1:1, j in -1:1]
    return Neighbourhood(hood)
end
function Base.getproperty(n::Neighbourhood, direction::Symbol)
    if direction == :south
        return getfield(n, :mask)[3,2]
    elseif direction == :south_east
        return getfield(n, :mask)[3,3]
    elseif direction == :east
        return getfield(n, :mask)[2,3]
    elseif direction == :north_east
        return getfield(n, :mask)[1,3]
    elseif direction == :north
        return getfield(n, :mask)[1,2]
    elseif direction == :north_west
        return getfield(n, :mask)[1,1]
    elseif direction == :west
        return getfield(n, :mask)[2,1]
    elseif direction == :south_west
        return getfield(n, :mask)[3,1]
    end
end
Base.size(::Neighbourhood) = (3,3)
Base.getindex(n::Neighbourhood, args...) = getindex(getfield(n, :mask), args...)

"""
    minimum_area_cycle_old(data, interval; threshold=Inf, step_limit=1_000_000)

NOTE: This is the original version of minimum_area_cycle, but runs slower.

Find the minimum area cycle in `data` for a given `interval`. The `data` and `threshold`
must be match the inputs to ripserer when computing persistent homology, and the `interval`
must be a zero-dimensional interval with a representative.  `step_limit` is a limit to the
length of the cycle returned and is used to prevent infinite loops.
"""
function minimum_area_cycle_old(interval; threshold=Inf, step_limit=1_000_000)

    # Mask is an array double the size of the image and mask[i,j]=true if (i/2,j/2) is in
    # the representative. We want to walk the boundary of the region that contains `true` in
    # counterclockwise order.
    orig_mask = representative_mask(interval, threshold)
    mask = zeros(Bool, size(orig_mask) .* 2)
    for i in axes(orig_mask, 1), j in axes(orig_mask, 2)
        mask[2i-1:2i, 2j-1:2j] .= orig_mask[i, j]
    end

    curr_index = CartesianIndex(findfirst(mask)) + CartesianIndex(-1, 0)
    start_index = curr_index

    boundary = Tuple{Float64, Float64}[]

    for _ in 1:step_limit
        neighbourhood = Neighbourhood(mask, curr_index)
        if length(boundary) > 3 && Tuple(curr_index) == boundary[end-1]
            @error "Cycle detected in $interval" neighbourhood
            break
        end
        @assert !get(mask, curr_index, false)
        point = fld.(Tuple(curr_index), 2) .+ 0.5
        # each point is visited twice, so we prevent adding duplicates here
        if get(boundary, lastindex(boundary), nothing) ≠ point
            push!(boundary, point)
        end
        # move S
        if (neighbourhood.west || neighbourhood.south_west) && !neighbourhood.south
            curr_index += CartesianIndex(1, 0)
        # move E
        elseif (neighbourhood.south || neighbourhood.south_east) && !neighbourhood.east
            curr_index += CartesianIndex(0, 1)
        # move N
        elseif (neighbourhood.east || neighbourhood.north_east) && !neighbourhood.north
            curr_index += CartesianIndex(-1, 0)
        # move W
        elseif (neighbourhood.north || neighbourhood.north_west) && !neighbourhood.west
            curr_index += CartesianIndex(0, -1)
        else
            # this should never be reached
            @assert false
        end

        if curr_index == start_index
            break
        end
    end
    return boundary
end

"""
    minimum_area_cycle(interval; threshold::Float64=Inf, step_limit::Int=1_000_000)
    
NOTE: This is a revised version of minimum_area_cycle_old that runs faster.

Find the minimum area cycle in `data` for a given `interval`. The `data` and `threshold`
must be match the inputs to ripserer when computing persistent homology, and the `interval`
must be a zero-dimensional interval with a representative.  `step_limit` is a limit to the
length of the cycle returned and is used to prevent infinite loops.

# Arguments
- `interval`: A PersistenceInterval with representative and birth_simplex
- `threshold`: Threshold value
- `step_limit`: Maximum number of steps in boundary tracing

# Returns
- Vector of (x, y) tuples at half-integer coordinates, 8-connected boundary (closed cycle)
- Set of CartesianIndex corresponding to pixels in the cycle
"""
function minimum_area_cycle(interval; threshold::Float64=Inf, step_limit::Int64=1_000_000)
    representative = interval.representative
    birth_simplex = interval.birth_simplex
    
    if isempty(representative)
        return Tuple{Float64,Float64}[]
    end
    
    # Build a map from CartesianIndex to cubes for efficient lookup
    pixel_to_cube = Dict(only(vertices(cube)) => cube for cube in representative)
    
    # Get starting pixel (global minimum - the birth simplex)
    birth_vertex = only(vertices(birth_simplex))
    
    # DFS to find connected component of sublevel set containing birth_vertex
    region_set = Set{CartesianIndex{2}}()
    stack = CartesianIndex{2}[birth_vertex]
    # allocate memory once
    sizehint!(region_set, length(representative)) 
    sizehint!(stack, 4length(representative))
    
    while !isempty(stack)
        curr = pop!(stack)
        
        curr in region_set && continue
        
        # Check if curr is in representative and has f < threshold
        cube = get(pixel_to_cube, curr, nothing)
        if isnothing(cube) || birth(cube) >= threshold
            continue
        end
        
        push!(region_set, curr)
        
        # Add 4-connected neighbours to stack
        for offset in OFFSETS
            neighbour = curr + offset
            neighbour ∉ region_set && push!(stack, neighbour)
        end
    end
    
    if isempty(region_set)
        return Tuple{Float64,Float64}[]
    end
    
    # Find starting point for tracing: leftmost-topmost pixel (lexicographically smallest)
    start_pixel_trace = minimum(x -> x.I, region_set)
    
    # Start one position west (outside region), on the 2x upscaled grid
    # Use plain integers instead of tuples
    curr_i = 2 * start_pixel_trace[1] - 1
    curr_j = 2 * start_pixel_trace[2] - 2
    start_i = curr_i
    start_j = curr_j
    
    boundary = Tuple{Float64,Float64}[]
    
    # Helper to check if a position on 2x grid is inside region
    # Inline everything to avoid allocations
    @inline function in_region_inline(i::Int, j::Int)
        pixel_idx = CartesianIndex(div(i + 1, 2), div(j + 1, 2))
        return pixel_idx in region_set
    end
    
    prev_point_i = typemin(Int)
    prev_point_j = typemin(Int)
    
    for _ in 1:step_limit
        # Get 8-neighbourhood - compute directly without allocations
        south = in_region_inline(curr_i + 1, curr_j)
        south_east = in_region_inline(curr_i + 1, curr_j + 1)
        east = in_region_inline(curr_i, curr_j + 1)
        north_east = in_region_inline(curr_i - 1, curr_j + 1)
        north = in_region_inline(curr_i - 1, curr_j)
        north_west = in_region_inline(curr_i - 1, curr_j - 1)
        west = in_region_inline(curr_i, curr_j - 1)
        south_west = in_region_inline(curr_i + 1, curr_j - 1)
        
        # Convert to half-integer coordinates
        # Use fld (floor division) to match minimum_area_cycle exactly
        point_i = fld(curr_i, 2)
        point_j = fld(curr_j, 2)
        
        # Avoid duplicates - only check coordinates directly
        if point_i != prev_point_i || point_j != prev_point_j
            push!(boundary, (point_i + 0.5, point_j + 0.5))
            prev_point_i = point_i
            prev_point_j = point_j
        end
        
        # Moore neighbourhood boundary tracing (same rules as minimum_area_cycle)
        # Move south
        if (west || south_west) && !south
            curr_i += 1
        # Move east
        elseif (south || south_east) && !east
            curr_j += 1
        # Move north
        elseif (east || north_east) && !north
            curr_i -= 1
        # Move west
        elseif (north || north_west) && !west
            curr_j -= 1
        else
            # Should not reach here for valid regions
            break
        end
        
        if curr_i == start_i && curr_j == start_j
            # Add the starting point again to close the cycle
            point_i = fld(curr_i, 2)
            point_j = fld(curr_j, 2)
            if point_i != prev_point_i || point_j != prev_point_j
                push!(boundary, (point_i + 0.5, point_j + 0.5))
            end
            break
        end
    end
    
    return (boundary, region_set)
end

"""
    shortest_and_minimum_area_cycles(data, interval; threshold=Inf)

Find the minimum area cycle and the shortest cycle in `data` for a given `interval`. The
`data` and `threshold` must be match the inputs to ripserer when computing persistent
homology, and the `interval` must be a zero-dimensional interval with a representative.
`step_limit` is a limit to the length of the cycle returned and is used to prevent infinite
loops.
"""
function shortest_and_minimum_area_cycles(data, interval; threshold=Inf)
    points = unique!(minimum_area_cycle(interval; threshold)[1])
    push!(points, points[1])
    tri = triangulate(unique(points))
    return points, points[get_convex_hull_vertices(tri)], area
end

function cycles_df(file::String; threshold, dir="../data")
    data = load_data(file, dir)
    return cycles_df(data, file; threshold)
end

function cycles_df(data, file=""; threshold)
    df = DataFrame()
    summary_df = DataFrame()

    neg_data = -data
    intervals = ripserer(
        Cubical(neg_data; threshold); dim_max=0, reps=true, verbose=true
    )[1]

    for interval in intervals
        area = length(interval.representative)
        int_tuple = (birth(interval), death(interval))
        max_pos = Tuple(only(vertices(interval.birth_simplex)))
        pers = persistence(interval)

        push!(summary_df, (; file, threshold, interval=int_tuple, persistence=pers, max_pos, area))

        cycle = minimum_area_cycle(interval; threshold)[1]

        for p in cycle
            push!(df, (
                ; file, threshold, interval=int_tuple, persistence=pers, max_pos, area, cycle=p
            ))
        end
    end
    return df, summary_df
end

function process_all(dir="../data", thresholds=-(0.02:0.02:0.1))
    for file in readdir(dir)
        cycle_df = DataFrame()
        summary_df = DataFrame()
        endswith(file, ".tif") || continue
        for threshold in thresholds
            @info "Proccing $file" threshold
            cyc, summ = cycles_df(file; dir, threshold)

            append!(cycle_df, cyc)
            append!(summary_df, summ)
        end
        f, _ = splitext(file)
        Arrow.write(joinpath(dir, "$f.cycles.arrow"), cycle_df; compress=:zstd)
        Arrow.write(joinpath(dir, "$f.summary.arrow"), summary_df; compress=:zstd)
    end

end
