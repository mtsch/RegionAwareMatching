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

# Example plots
if false
    data = load_data(2010)

    diag = ripserer(Cubical(-data); dim_max=0, merge_tree=true, verbose=true, reps=true)[1]
    diag2 = deepcopy(diag) # backup

    min_area = 2000
    filter_merge_tree!(x -> area(x) ≥ min_area, diag)
    filter_merge_tree!(x -> birth(x) < -0.1, diag)

    plot_cycles_merge_tree(data, diag)

    plot_all_at_threshold(data, diag, -0.5)
end
