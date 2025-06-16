include("data-loading.jl")
include("merge-trees.jl")
include("plotting-utils.jl")

function plot_variants(year)
    data = load_data(year)
    diag = ripserer(Cubical(-data; threshold=-0.05); dim_max=0, merge_tree=true, verbose=true, reps=true)[1]
    diag2 = deepcopy(diag)
    # Filter by persistence ≥ 0.05
    begin
        println("plotting a")
        diag = deepcopy(diag2)
        largest = argmax(area, diag)
        filtered = merge_subtree(diag, largest)

        filter_merge_tree!(x -> persistence(x) ≥ 0.05, filtered)
        @show length(filtered)
        @assert merge_tree_valid(filtered)

        f = plot_merge_tree_leaf_segmentation(
            data, filtered; title=L"%$year MTLS, threshold$=0.05$, persistence $\ge 0.05$"
        )
        save("../plots/$(year)_MTLS_pers.png", f)
    end

    # Filter by area > 2000
    begin
        println("plotting b")
        diag = deepcopy(diag2)
        largest = argmax(area, diag)
        filtered = merge_subtree(diag, largest)
        min_area = 2000
        filter_merge_tree!(x -> area(x) ≥ min_area, filtered)
        @show length(filtered)
        @assert merge_tree_valid(filtered)

        f = plot_merge_tree_leaf_segmentation(
            data, filtered; title=L"%$year MTLS, threshold$=0.05$, area $\ge 2000$"
        )
        save("../plots/$(year)_MTLS_area.png", f)
    end

    # No filter
#=
    begin
        println("plotting c")
        diag = deepcopy(diag2)
        largest = argmax(area, diag)
        filtered = merge_subtree(diag, largest)
        @show length(filtered)
        @assert merge_tree_valid(filtered)

        f = plot_merge_tree_leaf_segmentation(
            data, filtered; title=L"%$year MTLS, threshold$=0.05$"
        )
        save("../plots/$(year)_MTLS_peak.png", f)
    end
=#
end

if false
    for year in 1990:2020
        @show year
        plot_variants(year)
    end
end
