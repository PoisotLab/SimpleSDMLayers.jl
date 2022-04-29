# # Multivariate maps

using SimpleSDMLayers
using StatsPlots
using StatsPlots.PlotMeasures
using Statistics

# **Justification for this use case:** We can show more than one (specifically,
# two or three) variables on a single map, using a bivariate or trivariate color
# scale. In order to illustrate these mappings, we will look at the joint
# distribution of three measures: eveness of land use, terrain roughness, and
# proportion of urbanized land.

boundaries = (left=-12.0, right=30.0, bottom=36.0, top=72.0)
layer1 = convert(
    Float16,
    SimpleSDMPredictor(EarthEnv, HabitatHeterogeneity, 2; resolution=5, boundaries...),
)
layer2 = convert(
    Float16, SimpleSDMPredictor(EarthEnv, Topography, 7; resolution=5, boundaries...)
)

# The landcover data are finer than the other layers, so we will coarsen them.
# Because of rounding issues, the left and right cordinates need to be rounded.

layer3 = coarsen(
    convert(Float16, SimpleSDMPredictor(EarthEnv, LandCover, 9; boundaries...)),
    mean,
    (5, 5),
)
layer3.left = round(layer3.left; digits=0)
layer3.right = round(layer3.right; digits=0)

# We finally mask everything according to the first layer

layer2 = mask(layer1, layer2);
layer3 = mask(layer1, layer3);

# Note that bivariate maps usually work best when used with 9 classes in total
# (so 3 for each side). The next decision is to take a bivaraite color palette,
# and the combinations below are [commonly
# used](https://www.joshuastevens.net/cartography/make-a-bivariate-choropleth-map/).
# Note that you can definitely use [diverging
# colors](https://www.personal.psu.edu/cab38/ColorSch/Schemes.html) if you want.
# If you use colors in the RGBA format (*e.g.* `colorant"#ef0ce8c4"`), the color
# map will account for transparency. The color maps are actually passed as
# `ColorScheme` objects from the `ColorSchemes` package. The bivariate map
# itself is a call to plot. There are a number of color schemes defined in
# `SimpleSDMLayers.bivariates`. Internally, this will transform the layers into
# quantiles (determined by the `classes` keyword, defaults to 3):

plot(layer1, layer3; st=:bivariate, SimpleSDMLayers.bivariates.purple_yellow...)

# Note that you can use the `bivariate` shorthand as well:

pl1 = bivariate(layer1, layer3; classes=3, frame=:box, SimpleSDMLayers.bivariates.green_blue...)
xaxis!(pl1, "Longitude")
yaxis!(pl1, "Latitude")

# Note that *any* colorscheme is fair game, in case you want to mix diverging
# with directional:

bivariate(layer1, layer3; classes=12, frame=:box, grad1=ColorSchemes.tableau_orange_light, grad2=ColorSchemes.tableau_blue_green)
xaxis!("Longitude")
yaxis!("Latitude")

# We can repeat essentially the same process for the legend. In fact, inspecting
# the legend with the correct number of classes is probably a good idea in order
# to pick color schemes that mean something.

pl2 = bivariatelegend(layer1, layer3; classes=12, grad1=ColorSchemes.tableau_orange_light, grad2=ColorSchemes.tableau_blue_green)
xaxis!(pl2, layernames(EarthEnv, HabitatHeterogeneity, 2))
yaxis!(pl2, layernames(EarthEnv, LandCover, 9))

# And now, we can plot the legend next to the map - future releases of the
# package will hopefully offer this in a far more user friendly way.

plot(pl1, pl2; layout=@layout [a{0.75w} b])

# Using the `subplot` and `inset` arguments of Plots.jl, we can have the legend
# within the figure. Note how in this example we expand the limits on the x axis
# to make the legend fit, but also use more classes in the map to have a
# smoother result.

p1 = bivariate(
    layer1,
    layer3;
    classes=6,
    frame=:box,
    xlim=(-24, maximum(longitudes(layer1)))
)
xaxis!(p1, "Longitude")
yaxis!(p1, "Latitude")
p2 = bivariatelegend!(
    layer1,
    layer3;
    inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :left)),
    subplot=2,
    xlab=layernames(EarthEnv, HabitatHeterogeneity, 2),
    ylab=layernames(EarthEnv, LandCover, 9),
    guidefontsize=7
)

# Using a trivariate mapping follows the same process, with layers representing
# the red, green, and blue channel respectively.

plot(layer1, layer2, layer3; st=:trivariate, frame=:grid)

# There are two options for this type of plots. The first is `quantiles=true`
# (which maps quantiles rather than raw values), and the second is
# `simplex=false`, which makes all values sum to 1 within a pixel. For example:

trivariate(layer1, layer2, layer3; quantiles=true, simplex=true, frame=:grid)

# It is a good idea to question whether using `simplex` is appropriate. The
# legend can also be plotted using `trivariatelegend`:

trivariatelegend(layer1, layer2, layer3; quantiles=true, simplex=true)

# The legend function admits three additional arguments for the names of the
# `red`, `green`, and `blue` channels:

trivariatelegend(
    layer1,
    layer2,
    layer3;
    quantiles=true,
    simplex=true,
    red="Heterogeneous",
    green="Rough",
    blue="Urbanized"
)

# We can also combine the two elements:

tri1 = trivariate(
    layer1,
    layer2,
    layer3;
    xlim=(-24, maximum(longitudes(layer1))),
    frame=:grid,
    grid=false,
    simplex=false
)
xaxis!(tri1, "Longitude")
yaxis!(tri1, "Latitude")
tri2 = trivariatelegend!(
    layer1,
    layer2,
    layer3;
    inset=(1, bbox(0.01, 0.03, 0.35, 0.35, :top, :left)),
    subplot=2,
    red="Heterogeneity",
    green="Roughness",
    blue="Urban",
    annotationfontsize=6,
    simplex=false
)
