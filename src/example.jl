include("data-loading.jl")
include("merge-trees.jl")
include("plotting-utils.jl")

data = load_data(2010)

diag = ripserer(Cubical(-data); dim_max=0, merge_tree=true, verbose=true, reps=true)[1]

min_area = 2000
filtered = filter_merge_tree(x -> area(x) ≥ min_area, diag)
filter_merge_tree!(x -> birth(x) < -0.1, filtered)

plot_merge_tree_leaf_segmentation(data, filtered)

plot_all_at_threshold(data, filtered, -0.5)
