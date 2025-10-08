nclude("data-loading.jl")
include("merge-trees.jl")
include("plotting-utils.jl")

data = load_data(2010)

diag = ripserer(Cubical(-data); dim_max=0, merge_tree=true, verbose=true, reps=true)[1]
diag2 = deepcopy(diag)

# Filter by persistence ≥ 0.05
begin
    diag = deepcopy(diag2)
    largest = argmax(area, diag)
    filtered = merge_subtree(diag, largest)

    filter_merge_tree!(x -> persistence(x) ≥ 0.05, filtered)

    plot_merge_tree_leaf_segmentation(data, filtered; title=L"$year MTLS, persistence$≥ 0.05$")
end

# Filter by area > 2000 and birth > -0.05
begin
    diag = deepcopy(diag2)
    largest = argmax(area, diag)
    filtered = merge_subtree(diag, largest)
    min_area = 2000
    filter_merge_tree!(x -> area(x) ≥ min_area, filtered)
    filter_merge_tree!(x -> birth(x) ≤ -0.05, filtered)

    plot_merge_tree_leaf_segmentation(data, filtered; title=L"$year MTLS, peak$≥0.05$, area$≥2000$")
end

# Filter by birth > -0.05
begin
    diag = deepcopy(diag2)
    largest = argmax(area, diag)
    filtered = merge_subtree(diag, largest)
    filter_merge_tree!(x -> birth(x) ≤ -0.05, filtered)

    plot_merge_tree_leaf_segmentation(data, filtered; title=L"$year MTLS, peak$≥0.05$")
end

plot_all_at_threshold(data, filtered, -0.5)
