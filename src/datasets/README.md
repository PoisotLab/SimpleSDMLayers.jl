# Data provision interface

OK, this requires some notes, because adding new data sources has not been
enjoyable. Basically, the interface for a dataset is built around a
`LayerProvider` (*e.g.* WorldClim) and a `LayerDataset` (*e.g.* Elevation), the
combination of which defines a groups of rasters. The same `LayerDataset` can be
present within multiple `LayerProvider`s, and there is a function `provides` to
check that the dataset is provided by the layer. The first step in defining a
new provider is to overload the `provides` function to make it return `true` for
the relevant provider/dataset pairs. The second step is to define, for each
pair, the `layernames` method, which returns a tuple of strings describing what
each layer stores.

1. The "reference" implementation of the interface is the `WorldClim` provider,
   as it encompasses most of the features required to make it work

