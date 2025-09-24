# =============================================== #
# Utility functions for working with merge trees. #
# =============================================== #
using Ripserer, PersistenceDiagrams
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

Like `Base.filter!` but keeps merge tree structure intact. In-place version.

!!! warning
    Doesn't work if the parent of a kept interval is removed.
"""
function filter_merge_tree!(f, diagram)
    sort!(diagram, by=birth)

    curr_i = 0

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

Like `Base.filter` but keeps merge tree structure intact.

!!! warning
    Doesn't work if the parent of a kept interval is removed.
"""
function filter_merge_tree(f, diagram)
    result = deepcopy(diagram)
    return filter_merge_tree!(f, result)
end

"""
    merge_subtree(diagram, root)

Create a new diagram that includes the interval `root` and all its descendants.
"""
function merge_subtree(diagram, root)
    intervals = PersistenceInterval[]
    _merge_subtree!(intervals, root)
    return PersistenceDiagram(intervals, diagram.meta)
end

function _merge_subtree!(intervals, interval)
    push!(intervals, interval)
    for child in interval.children
        merge_subtree!(intervals, child)
    end
end
