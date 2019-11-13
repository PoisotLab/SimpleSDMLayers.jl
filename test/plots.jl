module SSLTestPlots
using SimpleSDMLayers
using Test
using Plots

temperature, precipitation = worldclim([1,12])

ispath("gallery") || mkpath("gallery")

heatmap(temperature, c=:RdYlBu_r, title="Temperature", frame=:box, clim=(-50,50))
xaxis!("Longitude")
yaxis!("Latitude")
savefig(joinpath("gallery", "heatmap.png"))

heatmap(coarsen(temperature, minimum, (10,10)), c=:RdYlBu_r, title="Temperature", frame=:box, clim=(-50,50))
xaxis!("Longitude")
yaxis!("Latitude")
savefig(joinpath("gallery", "heatmap_scaledown.png"))

histogram(precipitation, leg=false)
xaxis!("Precipitation")
savefig(joinpath("gallery", "density.png"))

plot(temperature, precipitation, leg=false, c=:grey, msc=:grey, alpha=0.5)
xaxis!("Temperature")
yaxis!("Precipitation")
savefig(joinpath("gallery", "scatter.png"))

end
