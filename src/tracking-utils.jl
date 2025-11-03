include("merge-trees.jl");
include("shortest-rep.jl");

using PersistenceDiagrams

# Track intervals

struct Vine{T<:Real}
	times::Vector{T}
	itvs::Vector{PersistenceInterval}
end
Base.length(vine::Vine) = length(vine.itvs)

function build_vines(years, vtx2itv_dicts, matching_prefix; info=false)
    # lookup uncompleted vines by most recent region, identified by birth vertex
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


# Track cycles

Cycle = Vector{Tuple{Float64,Float64}}

struct RegionVine{T<:Real}
	times::Vector{T}
	itvs::Vector{PersistenceInterval}
    cycles::Vector{Cycle}
    ischild_vec::Vector{Bool}
end
Base.length(vine::RegionVine) = length(vine.itvs)

# we only match two regions if they contain each other's minimum
# returns 3 booleans and a tuple:
# 1. whether child cycle of itv1 is matched
# 2. whether parent cycle of itv1 is matched
# 3. whether child-parent in itv1 switches to parent-child in itv2
# 4. 4-tuple of Union{Missing,Cycle}: (child cycle of itv1, parent cycle of itv1, child cycle of itv2, parent cycle of itv2)
#    missing if parent is nothing
function match_which(itv1, itv2)    
    if !isfinite(itv1.death) && !isfinite(itv2.death)
        cyc = minimum_area_cycle(itv1)[1]
        return (true, false, false, (cyc, missing, cyc, missing))
    end
    @assert isfinite(itv1.death) && isfinite(itv2.death)   

    child1_cyc, child1_reg = minimum_area_cycle(itv1)
    parent1_cyc, parent1_reg = minimum_area_cycle(itv1.parent; threshold=itv1.death)
    child2_cyc, child2_reg = minimum_area_cycle(itv2)
    parent2_cyc, parent2_reg = minimum_area_cycle(itv2.parent; threshold=itv2.death)
    cycles = (child1_cyc, parent1_cyc, child2_cyc, parent2_cyc)

    child1_vertex, parent1_vertex = only(vertices(itv1.birth_simplex)), only(vertices(itv1.parent.birth_simplex))
    child2_vertex, parent2_vertex = only(vertices(itv2.birth_simplex)), only(vertices(itv2.parent.birth_simplex))

    child1_child2 = (child1_vertex ∈ child2_reg) && (child2_vertex ∈ child1_reg)
    parent1_parent2 = (parent1_vertex ∈ parent2_reg) && (parent2_vertex ∈ parent1_reg)
    child1_parent2 = (child1_vertex ∈ parent2_reg) && (parent2_vertex ∈ child1_reg)
    parent1_child2 = (parent1_vertex ∈ child2_reg) && (child2_vertex ∈ parent1_reg)
    
    if (child1_child2 + parent1_parent2) >= (child1_parent2 + parent1_child2)
        return (child1_child2, parent1_parent2, false, cycles)
    else
        return (child1_parent2, parent1_child2, true, cycles)
    end
end

function build_region_vines(years, vtx2itv_dicts, matching_prefix; info=false)
    # lookup uncompleted vines by most recent region, identified by birth vertex and ischild
    active = Dict{Tuple{Int64,Int64,Bool}, RegionVine}(); 
    completed = RegionVine[]; # completed vines

    for yr1 in years[1:end-1]
        yr2 = yr1 + 1
        to_extend = Dict{Tuple{Int64,Int64,Bool}, RegionVine}(); # vines that can continue to be to_extend
        matching_yr1 = DataFrame(
            Arrow.Table(read("$(matching_prefix)_$(yr1)_$(yr2).arrow")),
        )
        for (left, right) in eachrow(matching_yr1)
            (ismissing(left) || ismissing(right)) && continue
            itv1 = vtx2itv_dicts[yr1][left]
            itv2 = vtx2itv_dicts[yr2][right]
            match_child1, match_parent1, switch, cycles = match_which(itv1, itv2)
            child1_cyc, parent1_cyc, child2_cyc, parent2_cyc = cycles
            GC.gc()

            for left_ischild in [true, false]
                !(left_ischild ? match_child1 : match_parent1) && continue # check if there is a reigon match
                right_ischild = switch ? !left_ischild : left_ischild
                key1 = (left..., left_ischild)
                key2 = (right..., right_ischild)

                vine = pop!(active, key1, RegionVine([yr1], [vtx2itv_dicts[yr1][left]], [left_ischild ? child1_cyc : parent1_cyc], [left_ischild]))
                push!(vine.times, yr2)
                push!(vine.itvs, vtx2itv_dicts[yr2][right])
                push!(vine.cycles, right_ischild ? child2_cyc : parent2_cyc)
                push!(vine.ischild_vec, right_ischild)
                to_extend[key2] = vine
            end        
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

### OLD VERSION

# maps interval to set of local minima contained inside representative
# function build_vertices_dict(diagram)
#     vertices_dict = Dict{PersistenceInterval, Set{CartesianIndex{2}}}()
#     for itv in sort(diagram.intervals, by=death)
#         vs = Set(vertices(itv.birth_simplex))
#         for child in filter(in(itv.children), diagram.intervals)
#             union!(vs, vertices_dict[child])
#         end
#         vertices_dict[itv] = vs
#     end
#     return vertices_dict
# end

# we only match two regions if they contain each other's local minima
# returns 3 booleans and a tuple:
# 1. whether child cycle of itv1 is matched
# 2. whether parent cycle of itv1 is matched
# 3. whether child-parent in itv1 switches to parent-child in itv2
# 4. 4-tuple of Union{Missing,Cycle}: (child cycle of itv1, parent cycle of itv1, child cycle of itv2, parent cycle of itv2)
#    missing if parent is nothing
# function match_which(itv1, itv2, yr1, yr2, itv2vertices_dicts)    
#     if !isfinite(itv1.death) && !isfinite(itv2.death)
#         cyc = minimum_area_cycle(itv1)[1]
#         return (true, false, false, (cyc, missing, cyc, missing))
#     end
#     @assert isfinite(itv1.death) && isfinite(itv2.death)   

#     child1_cyc, child1_reg = minimum_area_cycle(itv1)
#     parent1_cyc, parent1_reg = minimum_area_cycle(itv1.parent; threshold=itv1.death)
#     child2_cyc, child2_reg = minimum_area_cycle(itv2)
#     parent2_cyc, parent2_reg = minimum_area_cycle(itv2.parent; threshold=itv2.death)
#     cycles = (child1_cyc, parent1_cyc, child2_cyc, parent2_cyc)

#     vertices_dict1, vertices_dict2 = itv2vertices_dicts[yr1], itv2vertices_dicts[yr2]
#     child1_vertices, parent1_vertices = vertices_dict1[itv1], filter(in(parent1_reg), vertices_dict1[itv1.parent])
#     child2_vertices, parent2_vertices = vertices_dict2[itv2], filter(in(parent2_reg), vertices_dict2[itv2.parent])

#     child1_child2 = (child1_vertices ⊆ child2_reg) && (child2_vertices ⊆ child1_reg)
#     parent1_parent2 = (parent1_vertices ⊆ parent2_reg) && (parent2_vertices ⊆ parent1_reg)
#     child1_parent2 = (child1_vertices ⊆ parent2_reg) && (parent2_vertices ⊆ child1_reg)
#     parent1_child2 = (parent1_vertices ⊆ child2_reg) && (child2_vertices ⊆ parent1_reg)
    
#     if (child1_child2 + parent1_parent2) >= (child1_parent2 + parent1_child2)
#         return (child1_child2, parent1_parent2, false, cycles)
#     else
#         return (child1_parent2, parent1_child2, true, cycles)
#     end
# end

# function build_region_vines(years, vtx2itv_dicts, itv2vertices_dicts, matching_prefix; info=false)
#     # lookup uncompleted vines by most recent region, identified by birth vertex and ischild
#     active = Dict{Tuple{Int64,Int64,Bool}, RegionVine}(); 
#     completed = RegionVine[]; # completed vines

#     for yr1 in years[1:end-1]
#         yr2 = yr1 + 1
#         to_extend = Dict{Tuple{Int64,Int64,Bool}, RegionVine}(); # vines that can continue to be to_extend
#         matching_yr1 = DataFrame(
#             Arrow.Table(read("$(matching_prefix)_$(yr1)_$(yr2).arrow")),
#         )
#         for (left, right) in eachrow(matching_yr1)
#             (ismissing(left) || ismissing(right)) && continue
#             itv1 = vtx2itv_dicts[yr1][left]
#             itv2 = vtx2itv_dicts[yr2][right]
#             match_child1, match_parent1, switch, cycles = match_which(itv1, itv2, yr1, yr2, itv2vertices_dicts)
#             child1_cyc, parent1_cyc, child2_cyc, parent2_cyc = cycles

#             for left_ischild in [true, false]
#                 !(left_ischild ? match_child1 : match_parent1) && continue # check if there is a reigon match
#                 right_ischild = switch ? !left_ischild : left_ischild
#                 key1 = (left..., left_ischild)
#                 key2 = (right..., right_ischild)

#                 vine = pop!(active, key1, RegionVine([yr1], [vtx2itv_dicts[yr1][left]], [left_ischild ? child1_cyc : parent1_cyc], [left_ischild]))
#                 push!(vine.times, yr2)
#                 push!(vine.itvs, vtx2itv_dicts[yr2][right])
#                 push!(vine.cycles, right_ischild ? child2_cyc : parent2_cyc)
#                 push!(vine.ischild_vec, right_ischild)
#                 to_extend[key2] = vine
#             end        
#         end
#         append!(completed, collect(values(active)));
#         if info
#             @info "$yr1-$yr2" length(to_extend) length(completed)
#         end
#         active = to_extend
#         GC.gc()
#     end
#     append!(completed, collect(values(active))); # complete all active vines
#     return completed
# end

# interval to all birth vertices contained by interval (all local minima within cycle at death)
# itv2vertices_dicts = Dict(year => build_vertices_dict(diagrams[year]) for year in all_years);