using NCDatasets

file = "/workspaces/Fractal-Dimensions-of-Clouds/Data/dales_ensemble_final_masks.nc"

ds = NCDataset(file)
println(ds)
close(ds)