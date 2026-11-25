using NCDatasets
using FractalDimensions, ClimateBase, ComplexityMeasures
using CairoMakie
using Distances: Euclidean
using Statistics: mean

file1 = raw"/workspaces/Fractal-Dimensions-of-Clouds/Data/rlut_12am_feb3.nc"
file2 = raw"/workspaces/Fractal-Dimensions-of-Clouds/Data/rlut_12pm_feb3.nc"

ds1 = NCDataset(file1)
ds2 = NCDataset(file2)

println("Midnight: ", size(ds1["rlut"]))
println("Midday: ", size(ds2["rlut"]))

close(ds1)
close(ds2)