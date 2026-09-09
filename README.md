1.  Create a directory named `rasters`. From that directory, download the median *dhps540E* rasters from the [`dhps_mapping`](https://github.com/jflegg/dhps_mapping/tree/main) repository, for example, with the following Bash commands:
```
for year in $(seq 1990 2020); do
  curl -L -o Resistance_median_${year}.flt \
       https://raw.githubusercontent.com/jflegg/dhps_mapping/main/540%20rasters/Resistance_median_${year}.flt
  curl -L -o Resistance_median_${year}.hdr \
       https://raw.githubusercontent.com/jflegg/dhps_mapping/main/540%20rasters/Resistance_median_${year}.hdr
done
```

2.  Run `src/flt-to-tif-conversion.R` to convert the rasters to TIF files, followed by `src/extract-main-component.R` to extract only the largest connected component of each map.

3.  Run `src/region-wasmatch.jl` to produce persistence diagrams for all maps and perform region-aware matching.

4.  Figures 1 and 2 of the main text are produced by the scripts `src/pipeline-figures.jl` and `src/tracking-results.jl`, respectively.