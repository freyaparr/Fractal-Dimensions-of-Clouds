# Loading the data set
using CairoMakie, FractalDimensions, ClimateBase, DrWatson

cd("/workspaces/Fractal-Dimensions-of-Clouds/Correlation Dimension")

include("cloud_mask.jl")

file="/workspaces/Fractal-Dimensions-of-Clouds/Data/botany_2D_export.nc"

cloud = ncread(file, "lwp")
mask=makeMask(cloud, 0.01; below=false)

# making a grid of points
indices= findall(mask)
rows=[Float64(i[1]) for i in indices]
cols=[Float64(i[2]) for i in indices]
X=StateSpaceSet(rows, cols)
# range of values to test
εs= 10 .^(1:0.05:2)
# correlation sum method from FractalDimensions
Cs=correlationsum(X,εs;w=0,q=2)

# remove values of 0 and 1 for logs 
valid =(Cs.>0) .& (Cs.<1)
εs=εs[valid]
Cs=Cs[valid]

#
Δ, Δ05, Δ95 = slopefit(log2.(εs), log2.(Cs), LinearRegression())
# plotting 
fig, ax, obj = lines(log2.(εs), log2.(Cs))
ax.title = "slope = "*string(Δ)
ax.xlabel = "log2(ε)"
ax.ylabel = "log2(C_q(ε))"
save("CS_graph.png", fig)