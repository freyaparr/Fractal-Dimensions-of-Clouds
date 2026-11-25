using CairoMakie
using Statistics

#Identify the ensemble members with valid scaling exponents at BOTH times
valid = isfinite.(exponents1) .& isfinite.(exponents2)

xvals = exponents1[valid] #midnight
yvals = exponents2[valid] #midday

#correlation calulcation using the statistics package
r = cor(xvals, yvals)

#need the minimum and maximum values to see how big the plot is and determin how far the y=x line goes
low = min(minimum(xvals), minimum(yvals))
high = max(maximum(xvals), maximum(yvals))

fig = Figure()
ax = Axis(fig[1,1], xlabel = "Midnight ζ₂", ylabel = "Midday ζ₂", title = "Midday vs Midnight Radiative Scaling")

scatter!(ax, xvals, yvals, color = :blue, label = "Ensemble members")
lines!(ax, [low, high], [low, high], color = :black, linestyle = :dash, label = "y = x")

axislegend(ax, position = :rb)

println("Correlation = ", r)

save("rs_correlation_plot.png", fig)