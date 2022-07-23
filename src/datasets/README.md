# Data provision interface

OK, this requires some notes, because adding new data sources has not been
enjoyable. So there it goes. The data provision interface.

> The "reference" implementation of the interface is the `WorldClim` provider,
as it encompasses most of the features required to make it work: multiple
datasets, one of which has climate change scenarios, and a zip-based data origin
that shows why download/storage is important.

The goal of the interface from a maintainer point of view is to simplify the
addition of new data sources, which are identified by a `LayerProvider` and a
`LayerDataset`. The goal from a user point of view is to be able to call
`SimpleSDMPredictor` and get these data without needing to know where they are
stored, in which format, what the bounds are, *etc*. The endgoal is to make the
*user-facing* syntax simple, and agnostic to the format of the data.

The core of the interface for a new dataset is built around a `LayerProvider`
(this represents, more or less, a website on which there are layers we want to
use, *e.g.* WorldClim) and a `LayerDataset` (one thing you would get after a
hopefully finite number of clicks on the website; *e.g.* elevation).

The *unique* combination of a provider and a dataset defines a group of rasters,
which can have multiple layers, or can have multiple time points, or climate
scenarios. Part of the complexity is that the same `LayerDataset` can be present
within multiple `LayerProvider`s -- for example, WorldClim and CHELSA both
export a BioClim dataset (albeit at different resolutions, and with different
climate change models).

The first step in defining a dataset provider is to overload the `provides`
function -- by default, for any pair of provider/dataset, it will return
`false`. For any pair that is provided, it need to return `true`.

The second step in defining a dataset provider is to define, for each pair, the
`layernames` method, which returns **a tuple of strings** describing what each
layer stores.
