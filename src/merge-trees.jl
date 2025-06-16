using Ripserer, PersistenceDiagrams
using GLMakie
using LaTeXStrings


area(interval) = length(interval.representative)

"""
    merge_tree_valid(diagram)

Testing utility - check that merge tree is valid
"""
function merge_tree_valid(diagram)
    for int in diagram
        if !isnothing(int.parent)
            if int ∉ int.parent.children
                return false
            end
            if int.parent ∉ diagram
                return false
            end
        end
        for child in int.children
            child.parent ≡ int || return int, child
            child in diagram
        end
    end
    return true
end

"""
    filter_merge_tree!(f, diagram)

Like `filter!` but keeps merge tree structure intact. In-place version.
TODO: doesn't work if the parent of a kept interval is removed.
"""
function filter_merge_tree!(f, diagram)
    sort!(diagram, by=birth)

    curr_i = 0
    for i in eachindex(diagram)
        int = diagram[i]
        !f(int) && continue
    end

    for i in eachindex(diagram)
        int = diagram[i]
        if !f(int)
            if !isnothing(int.parent)
                child_idx = findfirst(x -> x ≡ int, int.parent.children)
                deleteat!(int.parent.children, child_idx)
            end
        else
            curr_i += 1
            diagram.intervals[curr_i] = int
        end
    end
    resize!(diagram.intervals, curr_i)
    return diagram
end

"""
    filter_merge_tree(f, diagram)

Like `filter` but keeps merge tree structure intact.
TODO: doesn't work if the parent of a kept interval is removed.
"""
function filter_merge_tree(f, diagram)
    result = deepcopy(diagram)
    return filter_merge_tree!(f, result)
end

function cycle_at_first_merge(int)
    if !isempty(int.children)
        threshold = minimum(death, int.children)
    else
        threshold = death(int)
    end
    return minimum_area_cycle(int; threshold), threshold
end


function get_subtree!(intervals, interval)
    push!(intervals, interval)
    for child in interval.children
        get_subtree!(intervals, child)
    end
end

function merge_subtree(diagram, root)
    intervals = PersistenceInterval[]
    get_subtree!(intervals, root)
    return PersistenceDiagram(intervals, diagram.meta)
end
