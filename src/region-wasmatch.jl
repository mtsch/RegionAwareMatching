include("data-loading.jl");
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

using Hungarian, LinearAlgebra, PersistenceDiagrams, StatsBase
using PersistenceDiagrams: Matching, MatchingDistance
import PersistenceDiagrams: _distance, _adjacency_matrix



### Region-aware matching, similar to https://arxiv.org/abs/2510.16486

struct RegionWasserstein <: MatchingDistance
    p::Float64
    q::Float64
    death_cutoff::Float64

    RegionWasserstein(p=2, q=2; death_cutoff=Inf) = new(Float64(p), Float64(q), Float64(death_cutoff))
end

function (rwas::RegionWasserstein)(
    left::PersistenceDiagram, right::PersistenceDiagram,
    left_raster::Raster, right_raster::Raster,
    left_seg_dict::Dict, right_seg_dict::Dict; 
	matching=false, 
)
    if length(left) == 0 & length(right) == 0
        if matching
            return Matching(left, right, 0, Pair{Int,Int}[], false)
        else
            return 0.0
        end
    end

    if count(!isfinite, left) == count(!isfinite, right)
        adj = _adjacency_matrix(right, left, right_raster, left_raster, right_seg_dict, left_seg_dict, rwas);
        match = collect(i => j for (i, j) in enumerate(hungarian(adj)[1]))

        distance = sum(adj[i, j] for (i, j) in match)^(1 / rwas.p)

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

function check_enclosure(cidx::CartesianIndex{2}, itv::PersistenceInterval, seg_dict::SegDict)
    if cidx ∈ seg_dict[itv]
        return true
    end
    return any(check_enclosure(cidx, child, seg_dict) for child in itv.children)
end

function _distance(
	left_itv::PersistenceInterval, right_itv::PersistenceInterval,
    left_raster::Raster, right_raster::Raster,
    left_seg_dict::SegDict, right_seg_dict::SegDict, 
	q::Float64, death_cutoff::Float64=Inf; verbose=false
)

    left_death = min(left_itv.death, death_cutoff)
    right_death = min(right_itv.death, death_cutoff)

    left_idxs = left_seg_dict[left_itv]
    right_idxs = right_seg_dict[right_itv]

    # regions need to enclose each other's birth vertex to be matched
    if !(
        check_enclosure(only(vertices(left_itv.birth_simplex)), right_itv, right_seg_dict) &&
        check_enclosure(only(vertices(right_itv.birth_simplex)), left_itv, left_seg_dict)
    )
        return Inf 
    end

    # shared = collect(left_idxs ∩ right_idxs)
    # other = collect(symdiff(left_idxs, right_idxs))
    # shared_diffs = getindex.(Ref(left_raster), shared) .- getindex.(Ref(right_raster), shared)  
    # diffs = [
    #     shared_diffs .- mean(shared_diffs);
    #     getindex.(Ref(left_raster), other) .- getindex.(Ref(right_raster), other);
    #     left_death === Inf && right_death === Inf ? 0. : left_death - right_death
    # ]

    combined = collect(left_idxs ∪ right_idxs)  
    # diffs = [
    #     getindex.(Ref(left_raster), combined) .- getindex.(Ref(right_raster), combined);
    #     left_death === Inf && right_death === Inf ? 0. : left_death - right_death
    # ]

    combined_diffs = getindex.(Ref(left_raster), combined) .- getindex.(Ref(right_raster), combined)
    diffs = [
        combined_diffs .- mean(combined_diffs);
        left_itv.birth - right_itv.birth;
        left_death === Inf && right_death === Inf ? 0. : left_death - right_death
    ]

    if verbose
        display(lines(abs.(diffs)))
    end
    return norm(diffs, q)
end

function _diagonal_distance(
    itv::PersistenceInterval, raster::Raster, seg_dict::SegDict, q::Float64, death_cutoff::Float64=Inf
)
    # N.B. work in negated space
    vals = -getindex.(Ref(raster), seg_dict[itv])
    itv_death = min(itv.death, death_cutoff)
    mid_val = (itv.birth + itv_death)/2
    # mid_val = mean(vals)
    isfinite(itv_death) && push!(vals, itv.death) # See Eqn. (3) of https://arxiv.org/abs/2510.16486
    return norm(vals .- mid_val, q)
end

function _adjacency_matrix(
	left::PersistenceDiagram, right::PersistenceDiagram,
    left_raster::Raster, right_raster::Raster,
    left_seg_dict::SegDict, right_seg_dict::SegDict, 
	rwas::RegionWasserstein
)
    n = length(left)
    m = length(right)
    adj = fill(Inf, n + m, m + n)

    for j in eachindex(left), i in eachindex(right)
        adj[i, j] = _distance(
            left[j], right[i], left_raster, right_raster, 
            left_seg_dict, right_seg_dict, rwas.q, rwas.death_cutoff
        )
    end

    for i in 1:n
        adj[i + m, i] = _diagonal_distance(left[i], left_raster, left_seg_dict, rwas.q, rwas.death_cutoff)
    end
    for j in 1:m
        adj[j, j + n] = _diagonal_distance(right[j], right_raster, right_seg_dict, rwas.q, rwas.death_cutoff)
    end
    adj[(m + 1):(m + n), (n + 1):(n + m)] .= 0.0

    if rwas.p ≠ 1
        adj .^= rwas.p
    end
    return adj
end



## produces region-aware matchings for all pairs of consecutive years
function main()
    regenerate_diagrams() # comment out if diagrams already generated

    years = 1990:2020
    death_cutoff = -0.05
    pers_cutoff = 0.01
    cutoff_func(x) = persistence(x) >= pers_cutoff && birth(x) <= death_cutoff
    for yr1 in years[1:end-1]
        yr2 = yr1 + 1
        diagram1 = load_diagram(yr1; lazy=false); diagram2 = load_diagram(yr2; lazy=false);
        raster1 = load_data(yr1; raster=true, lazy=false); raster2 = load_data(yr2; raster=true, lazy=false);

        filter!(cutoff_func, diagram1.intervals)
        filter!(cutoff_func, diagram2.intervals)   
        for itv in diagram1.intervals
            filter!(cutoff_func, itv.children)
        end
        for itv in diagram2.intervals
            filter!(cutoff_func, itv.children)
        end     

        segmentation1 = build_segmentation(diagram1, raster1; death_cutoff);
        segmentation2 = build_segmentation(diagram2, raster2; death_cutoff);

        @info "$yr1-$yr2" length(diagram1) length(diagram2)

        region_matching = RegionWasserstein(2, 2; death_cutoff)(
            diagram1, diagram2, raster1, raster2, segmentation1, segmentation2; 
            matching=true
        );

        output = matching_to_df(region_matching)
        Arrow.write(joinpath(
            @__DIR__, "../data/match",
            "region_match_pers$(pers_cutoff)_death$(-death_cutoff)_$(yr1)_$(yr2).arrow"), output
        )
    end
    @info "fin."
end

if !isinteractive()
    main()
end
