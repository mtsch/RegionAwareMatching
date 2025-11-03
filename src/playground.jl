include("data-loading.jl");
include("merge-trees.jl");
include("plotting-utils.jl");
include("wasmatch.jl");

using StatsBase

all_years = 1990:2020

# file_postfix = "_cutoff0.0001"
file_postfix = "_area1000"

DATA_DIR = joinpath(@__DIR__, "../data/")
IMG_DIR = joinpath(@__DIR__, "../imgs/")
LWD = 1.5

avoid_colors = sequential_palette(240, 500)[[250, 500]]
palette = [
	RGB(0,0,0);
	distinguishable_colors(19, [RGB(0,0,0); RGB(1,1,1); avoid_colors], dropseed=true)
]

@time rasters = Dict(year => load_data(year; raster=true, lazy=false) for year in all_years);
@time diagrams = Dict(year => load_diagram(year; lazy=false) for year in all_years);

@time matchings = Dict(year => DataFrame(
	Arrow.Table(read("$DATA_DIR/match/match_$(year-1)_$(year)$(file_postfix).arrow")),
) for year in all_years[2:end]);

@time vtx2itv_dicts = Dict(year => Dict(
	birth_vertex(itv) => itv for itv in diagrams[year]
) for year in all_years);

struct Vine{T<:Real}
	times::Vector{T}
	itvs::Vector{PersistenceInterval}
end
Base.length(vine::Vine) = length(vine.itvs)

active = Dict{Tuple{Int64,Int64}, Vine}(); # lookup uncompleted vines by latest birth vertex
completed = Vine[];
for year in all_years[2:end]
	new_active = Dict{Tuple{Int64,Int64}, Vine}();
	for (left, right) in eachrow(matchings[year])
		ismissing(left) && continue
		if haskey(active, left) # vine already exists
			vine = active[left]
			if ismissing(right) # vine ends here
				push!(completed, vine)
			else # extend the vine
				push!(vine.times, year)
				push!(vine.itvs, vtx2itv_dicts[year][right])
				new_active[right] = vine
			end
		else # vine does not exist
			ismissing(right) && continue
			# create new vine
			new_active[right] = Vine([year-1, year], [vtx2itv_dicts[year-1][left], vtx2itv_dicts[year][right]])
		end
	end
	active = new_active
end
append!(completed, collect(values(active))); # complete all active vines
sort(countmap(length.(completed)))

itv2vine = Dict(
	(year, itv) => vine for vine in completed for (year, itv) in zip(vine.times, vine.itvs)
); # map (year, birth vertex) to vine
length(itv2vine)

# sorted_vines = sort(completed, by=(vine)->sum(persistence.(vine.itvs)), rev=true);
# vine = sorted_vines[3];
# persistence.(vine.itvs)
# mkpath(joinpath(IMG_DIR, "totpers3"));

sorted_vines = sort(completed, by=(vine)->sum((itv)->length(itv.representative), vine.itvs), rev=true);
vine = sorted_vines[3];
map((itv)->length(itv.representative), vine.itvs)
OUT_DIR = mkpath(joinpath(IMG_DIR, "area1000_totarea3"));

# for (year, itv) in zip(vine.times, vine.itvs)
# 	println("Original")
# 	@time c1 = minimum_area_cycle_old(itv.parent; threshold=itv.death)
# 	println("Rewrite")
# 	@time c2, _ = minimum_area_cycle(itv.parent; threshold=itv.death)
# 	@assert length(c1) == length(c2)
# 	@assert Set(c1) == Set(c2)
# end

for (year, itv) in zip(vine.times, vine.itvs)
	println(year)
	fig = Figure()
	ax = plot_heatmap!(fig, rasters[year]; title=string(year))
	plot_cycle!(ax, rasters[year], itv; birth_simplex=false, linewidth=LWD, color=palette[2])
	plot_cycle!(
		ax, rasters[year], itv.parent; 
		threshold=itv.death, linestyle=:dash, birth_simplex=false, linewidth=LWD, color=palette[2]
	)
	save(joinpath(OUT_DIR, "vine_$(year).png"), fig)
	# display(fig)
end

# plot top persistence cycles of each year
OUT_DIR = mkpath(joinpath(IMG_DIR, "top10cycles"))
for year in all_years
	println(year)
	fig = Figure()
	ax = plot_heatmap!(fig, rasters[year]; title=string(year))
	sorted_itvs = sort(diagrams[year],by=persistence, rev=true)
	for i in 1:10
		itv = sorted_itvs[i]
		plot_cycle!(ax, rasters[year], itv; birth_simplex=false, linewidth=LWD, color=palette[i])
	end
	save(joinpath(OUT_DIR, "$(year).png"), fig)
end

vine = sorted_vines[2]
# Investigate 1992-1993, 1997-1998, 1998-1999, 2008-2009, 2010-2011, 2014-2015, 2017-2018-2019, 2019-2020

# 2008-2009

# plot top persistence cycles of one year

year = 2008
sorted_itvs = sort(diagrams[year],by=persistence, rev=true);
sorted_itvs[5].parent === sorted_itvs[2]
sorted_itvs[2].parent === sorted_itvs[1]

year = 2009
sorted_itvs = sort(diagrams[year],by=persistence, rev=true);
sorted_itvs[3].parent === sorted_itvs[1]
sorted_itvs[2].parent === sorted_itvs[1]

# 2010-2011

year = 2010
sorted_itvs = sort(diagrams[year],by=persistence, rev=true);
sorted_itvs[3].parent === sorted_itvs[1]
sorted_itvs[2].parent === sorted_itvs[1]

year = 2011
sorted_itvs = sort(diagrams[year],by=persistence, rev=true);
sorted_itvs[7].parent === sorted_itvs[2]
sorted_itvs[2].parent === sorted_itvs[1]

# 2019-2020

vine.itvs[end-1]
vine.itvs[end]
sort(diagrams[2020],by=persistence, rev=true)[1:10]

birth_vertex(vine.itvs[end]) == birth_vertex(sort(diagrams[2020],by=persistence, rev=true)[8])

alt_itv = sort(diagrams[2020],by=persistence, rev=true)[2]
alt_vine = itv2vine[(2020, alt_itv)]

PersistenceDiagrams._distance(vine.itvs[end-1], vine.itvs[end])+PersistenceDiagrams._distance(alt_vine.itvs[end-1], alt_vine.itvs[end])

PersistenceDiagrams._distance(vine.itvs[end-1], alt_vine.itvs[end])+PersistenceDiagrams._distance(alt_vine.itvs[end-1], vine.itvs[end])

# follow main vine and alt vine
begin
	fig = Figure()
	year = 2019
	ax = plot_heatmap!(fig, rasters[year]; title="$year")
	sorted_itvs = sort(diagrams[year],by=persistence, rev=true)
	plot_cycle!(ax, rasters[year], vine.itvs[end-1]; birth_simplex=false, linewidth=LWD, color=palette[2], alpha=0.8)
	plot_cycle!(ax, rasters[year], alt_vine.itvs[end-1]; birth_simplex=false, linewidth=LWD, color=palette[3], alpha=0.8)
	fig
end

begin
	fig = Figure()
	year = 2020
	ax = plot_heatmap!(fig, rasters[year]; title="$year")
	sorted_itvs = sort(diagrams[year],by=persistence, rev=true)
	plot_cycle!(ax, rasters[year], vine.itvs[end]; birth_simplex=false, linewidth=LWD, color=palette[2], alpha=0.8)
	plot_cycle!(ax, rasters[year], alt_vine.itvs[end]; birth_simplex=false, linewidth=LWD, color=palette[3], alpha=0.8)
	fig
end

# plot top persistence cycles of one year
begin
	year = 2019
	fig = Figure()
	ax = plot_heatmap!(fig, rasters[year]; title=string(year))
	sorted_itvs = sort(diagrams[year],by=persistence, rev=true)
	for i in 2:8
		itv = sorted_itvs[i]
		plot_cycle!(ax, rasters[year], itv; birth_simplex=false, linewidth=LWD, color=palette[i])
	end
	fig
end

# 2019: 2 and 7?
year = 2019
sorted_itvs = sort(diagrams[year],by=persistence, rev=true);
sorted_itvs[1:10]
sorted_itvs[2].parent === sorted_itvs[1]
sorted_itvs[7].parent === sorted_itvs[2]

begin
	year = 2019
	fig = Figure()
	ax = plot_heatmap!(fig, rasters[year]; title=string(year))
	sorted_itvs = sort(diagrams[year],by=persistence, rev=true)
	i = 7
	itv = sorted_itvs[i]
	plot_cycle!(ax, rasters[year], itv; birth_simplex=false, linewidth=LWD, color=palette[i])
	plot_cycle!(ax, rasters[year], itv.parent; birth_simplex=false, linewidth=LWD, color=palette[i], threshold=itv.death, linestyle=:dash)

	fig
end

# 2020: 2 and 4?
year = 2020
sorted_itvs = sort(diagrams[year],by=persistence, rev=true);
sorted_itvs[1:10]
sorted_itvs[2].parent == sorted_itvs[1]
sorted_itvs[4].parent == sorted_itvs[1]

begin
	year = 2020
	fig = Figure()
	ax = plot_heatmap!(fig, rasters[year]; title=string(year))
	sorted_itvs = sort(diagrams[year],by=persistence, rev=true)
	i = 4
	itv = sorted_itvs[i]
	plot_cycle!(ax, rasters[year], itv; birth_simplex=false, linewidth=LWD, color=palette[i])
	plot_cycle!(ax, rasters[year], itv.parent; birth_simplex=false, linewidth=LWD, color=palette[i], threshold=itv.death, linestyle=:dash)

	fig
end

mkpath(joinpath(IMG_DIR, "area1000_first_merge_cycles"))
mkpath(joinpath(IMG_DIR, "area1000_all_cycles"))

for year in 1990:2020
	diagram = diagrams[year]
	itvs = filter((itv)->length(itv.representative)>=1000, diagram.intervals)
	itvs = sort(itvs, by=persistence, rev=true)

	# itvs = sort(diagram.intervals, by=persistence, rev=true)[1:20]
	# min_area = minimum(length ∘ representative, itvs)
	# println("$year: $min_area")
	# display(map(length ∘ representative, itvs))
	# break

	ncycles = 0

	fig = Figure()
	ax = plot_heatmap!(fig, rasters[year]; title=string(year))
	for i in 1:min(20, length(itvs))
		itv = itvs[i]
		children_in_itvs = filter((child)->child.parent == itv, itvs)
		# threshold = minimum(death, children_in_itvs; init=death(itv))
		# display((birth(itv), threshold))
		# plot_cycle!(
		# 	ax, rasters[year], itv; 
		# 	birth_simplex=false, threshold=threshold,
		# 	linewidth=LWD, color=palette[i]
		# )
		for child in children_in_itvs
			plot_cycle!(
				ax, rasters[year], itv; 
				birth_simplex=false, threshold=child.death,
				linewidth=LWD, color=palette[i]
			)
			ncycles += 1
		end
		plot_cycle!(
			ax, rasters[year], itv; 
			birth_simplex=false,
			linewidth=LWD, color=palette[i]
		)
		ncycles += 1		
	end
	println("$year: $(length(itvs)) $ncycles")
	display(fig)
	save(joinpath(IMG_DIR, "area1000_all_cycles/$(year).png"), fig)
	# break
end

# Old playground

getproperty.(completed, :times)
sort(countmap(length.(completed)))

longest = argmax(length, completed)
longest.itvs

birth_vertex(diagrams[1990][end])
row = findfirst(matchings[1991].left) do id
	!ismissing(id) && id == birth_vertex(diagrams[1990][end])
end
matchings[1991][1,:]

findfirst(completed) do vine
	vine.times[1] == 1990 && birth_vertex(vine.itvs[1]) == (1260, 1232)
end

length(completed[12989])
itv = diagrams[1990][1]
itv.parent

vine.itvs[end-1]
vine.itvs[end]

res1
