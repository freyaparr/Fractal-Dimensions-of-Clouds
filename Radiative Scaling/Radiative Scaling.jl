using NCDatasets
using FractalDimensions, ClimateBase, ComplexityMeasures
using Dates: hour
using CairoMakie
using Distances: Euclidean

# Botany data - Load the botany dataset using the file path on your computer
filename = raw"/workspaces/Fractal-Dimensions-of-Clouds/dales_radiation_rlut_data.nc"

file = filename
data = ncread(file, "rlut")
x, y = dims(data)

times = NCDataset(file)["time"][]
hours = hour(times)

euclidean = Euclidean()

# %%
# range of distances in km units
ks = 10 .^ (1:0.05:2) # equispaced exponentially
#ks = :0.5:150 # equispaced linearly
encoding = RectangularBinEncoding(FixedRectangularBinning(ks, 1, true))
δΦsum = zeros(length(ks) - 1)
counts = zeros(Int, length(ks) - 1)

function disc_sampling(x,y,R)
    # picks a random point
    i = rand(eachindex(x))
    j = rand(eachindex(y))

    # Generate a random distance and angle within the circle of radius R
    radius = R * 1000 * sqrt(rand()) 
    θ = 2π * rand() 

    # Calculate the new point's coordinates
    x2 = x[i] + radius * cos(θ)
    y2 = y[j] + radius * sin(θ)

    # Find the nearest grid indices
    i2 = argmin(abs.(x .- x2))
    j2 = argmin(abs.(y .- y2))
    return i, j, i2, j2
end

Nsample = 10000

for i in 5:100
    samples = [disc_sampling(x,y,i) for _ in 1:Nsample]
    for (i, j, i2, j2) in samples
        # distance in kilometers
        distance = euclidean(
        (x[i], y[j]),
        (x[i2], y[j2]) ) /1000

        bin_index = encode(encoding, distance)
        bin_index == -1 && continue # distance that doesn't fall into any bin

        Φ1 = data[i, j]
        Φ2 = data[i2, j2]
        δΦ = abs2(Φ1 - Φ2)

        ismissing(δΦ) && continue
        δΦsum[bin_index] += δΦ
        counts[bin_index] += 1
    end
end

# %% Plot and compute FD
using Statistics: mean

δΦs = δΦsum ./ counts
midpoints = (ks[1:end-1] + ks[2:end])/2

Δ, Δ05, Δ95 = slopefit(log10.(midpoints), log10.(δΦs), LinearRegression())

fig, ax, obj = lines(midpoints, δΦs)
lines!(ax, midpoints, minimum(δΦs) .* (midpoints .^ Δ); linestyle = :dash)
ax.xscale = log10
ax.yscale = log10
#xlims!(ax, 10, 10^1.9)

# also plot the colletced point density, to see how well our sampling strategy covers distance bins
ax2, obj = barplot(fig[2,1], midpoints, counts; width = diff(ks))
ax2.xscale = log10
ax2.yscale = log10
#xlims!(ax2, 10^1.1, 10^1.9)

ax.title = "slope = "*string(Δ)
ax2.xlabel = "r [km]"
ax.ylabel = "δΦ [W/m²]"
ax2.ylabel = "counts"

save("RadScal_graph.png", fig)
