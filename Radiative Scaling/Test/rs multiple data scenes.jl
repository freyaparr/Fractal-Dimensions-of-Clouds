using NCDatasets
using FractalDimensions, ClimateBase, ComplexityMeasures
using CairoMakie
using Distances: Euclidean
using Statistics: mean

# Botany data - Load the botany dataset using the file path on your computer
filename = raw"/workspaces/Fractal-Dimensions-of-Clouds/Data/dales_radiation_rlut_data_60hours.nc"

file = filename
data = ncread(file, "rlut")
x, y, times = dims(data)
#grid spacing (in m)
dx = x[2] - x[1]
dy = y[2] - y[1]
#total length of data in metres
#(we do know this from the report but this makes for better reproducibility if we use other data)
Lx = length(x) *dx
Ly = length(y) *dy

euclidean = Euclidean()

# %%
# range of distances in km units
ks = 10 .^ (0:0.1:2) # equispaced exponentially
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

    #wrap the points
    x2 = mod(x2, Lx)
    y2 = mod(y2, Ly) #this is simpler than using an if >Lx subtraction thing

    # Find the nearest grid indices
    i2 = argmin(abs.(x .- x2))
    j2 = argmin(abs.(y .- y2))
    return i, j, i2, j2
end

time_indices = 1:1:60
exponents = zeros(length(time_indices)) #this is to store the scaling expontents zeta_2 we get at each time frame
Nsample = 10000
#more efficient sampling so we dont generate new pairs for every member
all_samples = [[disc_sampling(x,y,i) for _ in 1:Nsample] for i in 1:150]

#basically chuck the code we had for the individual scene into a loop through each time (1 to 60)
for (k, t) in enumerate(time_indices)
    δΦsum .= 0
    counts .= 0

    for samples in all_samples
        for (i, j, i2, j2) in samples
            #new distance calculation that takes into accoiunt periodic data
            #(basically instead of 150km and 2km being 148km away, they're now only like 5km away)
            #picking which distance option is shortest 
            Δx = min(abs(x[i2]-x[i]), Lx-abs(x[i2]-x[i]))
            Δy = min(abs(y[j2]-y[j]), Ly-abs(y[j2]-y[j]))
            # distance in kilometers
            distance = euclidean((0, 0), (Δx, Δy)) / 1000

            bin_index = encode(encoding, distance)
            bin_index == -1 && continue # distance that doesn't fall into any bin

            Φ1 = data[i, j, t]
            Φ2 = data[i2, j2, t]
            δΦ = abs2(Φ1 - Φ2)

            ismissing(δΦ) && continue
            isnan(δΦ) && continue
            δΦsum[bin_index] += δΦ
            counts[bin_index] += 1
        end
    end


    # %% Plot and compute FD
    δΦs = δΦsum ./ counts
    midpoints = (ks[1:end-1] + ks[2:end])/2

    Δ, Δ05, Δ95 = slopefit(log10.(midpoints), log10.(δΦs), LinearRegression())
    exponents[k] = Δ
end

#replaced old plotting because we now want to jsut do a simple plot of the scaling exponents against times
hours = time_indices .- 1
fig=Figure()
ax = Axis(fig[1,1], xlabel = "Simulation Hour", ylabel = "Scaling Exponent ζ₂", title = "Radiative Scaling Over Time")

scatterlines!(ax, hours, exponents, color = :blue)

save("Radiative_Scaling_Over_Time.png", fig) 

#We have like 57.6 million sampled pairs so this will take agesss to run 
