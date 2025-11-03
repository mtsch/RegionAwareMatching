include("data-loading.jl")
include("merge-trees.jl");
include("shortest-rep.jl");

# copied from wasmatch.jl
function birth_vertex(interval)
    # the matching code will create new intervals on the diagonal to match to. This replaces
    # them with missing
    if hasproperty(interval, :birth_simplex)
        return to_indices(interval.birth_simplex)
    else
        return missing
    end
end

# copied from wasmatch.jl
function matching_to_df(m)
    result = DataFrame()
    for (left, right) in matching(m)
        push!(result, (left=birth_vertex(left), right=birth_vertex(right)); promote=true)
    end
    return result
end

using Hungarian, PersistenceDiagrams
using PersistenceDiagrams: Matching, MatchingDistance, _distance, _diagonal_interval
import PersistenceDiagrams: _distances, _adjacency_matrix

struct SpatialWasserstein <: MatchingDistance
    p::Float64
    q::Float64
    min_overlap::Float64

    SpatialWasserstein(p=1, q=Inf, min_overlap=0.0) = new(Float64(p), Float64(q), Float64(min_overlap))
end

RegionDict = Dict{Tuple{PersistenceInterval,Symbol},Set{CartesianIndex{2}}}

function build_region_dict(diagram::PersistenceDiagram)
    region_dict = Dict(
        (itv, :parent) => minimum_area_cycle(itv.parent; threshold=itv.death)[2]
        for itv in diagram if itv.parent !== nothing
    )
    merge!(region_dict, Dict(
        (itv, :child) => Set(map(only ∘ vertices, itv.representative))
        for itv in diagram if itv.parent !== nothing
    ))
    return region_dict
end

function jaccard_sim(x1, x2)
    smaller, larger = length(x1) <= length(x2) ? (x1, x2) : (x2, x1)
    intersection_size = count(in(larger), smaller)
    union_size = length(x1) + length(x2) - intersection_size
    return intersection_size / union_size
end

function _distances(
    left, right, 
    left_region_dict::RegionDict, right_region_dict::RegionDict, 
    p::Float64=1.0, q::Float64=Inf, min_overlap::Float64=0.0
)
    dists = zeros(length(right), length(left))
    for j in eachindex(left), i in eachindex(right)
        pers_dist = _distance(left[j], right[i], q)
        if !isfinite(left[j]) || !isfinite(right[i])
            dists[i, j] = pers_dist
        else
            left_child_reg = left_region_dict[(left[j], :child)]
            left_parent_reg = left_region_dict[(left[j], :parent)]
            right_child_reg = right_region_dict[(right[i], :child)]
            right_parent_reg = right_region_dict[(right[i], :parent)]
            # min captures worst overlap out of two pairs of regions
            # max captures best pairing out of two possible ways of pairing
            region_sim = max(
                min(jaccard_sim(left_child_reg, right_child_reg), jaccard_sim(left_parent_reg, right_parent_reg)), 
                min(jaccard_sim(left_child_reg, right_parent_reg), jaccard_sim(left_parent_reg, right_child_reg))
            )
            # cost if both intervals are matched to diagonal
            pers_diag_dist = (_distance(left[j], _diagonal_interval(left[j]), q)^p + _distance(right[i], _diagonal_interval(right[i]), q)^p)^(1/p)
            dists[i, j] = max(
                pers_dist,
                (pers_dist - pers_diag_dist) / (1 - min_overlap) * (region_sim - 1) + pers_dist
            )
        end
    end
    return dists
end

function _adjacency_matrix(
	left::PersistenceDiagram, right::PersistenceDiagram, 
	left_region_dict::RegionDict, right_region_dict::RegionDict, 
	power=1, q=Inf, min_overlap::Float64=0.0
)
    n = length(left)
    m = length(right)
    adj = fill(Inf, n + m, m + n)

    dists = _distances(left, right, left_region_dict, right_region_dict, power, q, min_overlap)
    adj[axes(dists)...] .= dists

    for i in 1:n
        adj[i + m, i] = _distance(left[i], _diagonal_interval(left[i]), q)
    end
    for j in 1:m
        adj[j, j + n] = _distance(right[j], _diagonal_interval(right[j]), q)
    end
    adj[(m + 1):(m + n), (n + 1):(n + m)] .= 0.0

    if power ≠ 1
        adj .^= power
    end
    return adj
end

function (w::SpatialWasserstein)(
    left::PersistenceDiagram, right::PersistenceDiagram; 
	matching=false, 
	left_region_dict::RegionDict=build_region_dict(left),
	right_region_dict::RegionDict=build_region_dict(right)
)
	if length(left) == 0 & length(right) == 0
        if matching
            return Matching(left, right, 0, Pair{Int,Int}[], false)
        else
            return 0.0
        end
    end

    if count(!isfinite, left) == count(!isfinite, right)
        adj = _adjacency_matrix(right, left, right_region_dict, left_region_dict, w.p, w.q, w.min_overlap)
        match = collect(i => j for (i, j) in enumerate(hungarian(adj)[1]))

        distance = sum(adj[i, j] for (i, j) in match)^(1 / w.p)

        if matching
            return Matching(left, right, distance, match, false)
        else
            return distance
        end
    else
        if matching
            return Matching(left, right, Inf, Pair{Int,Int}[], false)
        else
            return Inf
        end
    end
end

## produces matchings for all pairs of consecutive years
function main()
    years = 1990:2020
    for yr1 in years[1:end-1]
        yr2 = yr1 + 1
        diagram1 = load_diagram(yr1; lazy=false); diagram2 = load_diagram(yr2; lazy=false);

        # area_cutoff = 1000
        # filter!(x -> area(x) >= area_cutoff, diagram1.intervals)
        # filter!(x -> area(x) >= area_cutoff, diagram2.intervals)
        pers_cutoff = 0.01
        cutoff_func(x) = persistence(x) >= pers_cutoff
        filter!(cutoff_func, diagram1.intervals)
        filter!(cutoff_func, diagram2.intervals)   
        for itv in diagram1.intervals
            filter!(cutoff_func, itv.children)
        end
        for itv in diagram2.intervals
            filter!(cutoff_func, itv.children)
        end     

        @info "$yr1-$yr2" length(diagram1) length(diagram2)

        regdict1 = build_region_dict(diagram1);
        regdict2 = build_region_dict(diagram2);

        spatial_matching = SpatialWasserstein()(
            diagram1, diagram2; 
            left_region_dict=regdict1, right_region_dict=regdict2, matching=true
        );

        output = matching_to_df(spatial_matching)
        Arrow.write(joinpath(@__DIR__, "../data/match/spatial_match_pers$(pers_cutoff)_$(yr1)_$(yr2).arrow"), output)
    end
    @info "fin."
end

if !isinteractive()
    main()
end