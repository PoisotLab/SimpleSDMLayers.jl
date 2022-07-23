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
layer stores. This method is crucial because providers usually refer to layers
as codes or integer, and we want users to be able to figure out which layer is
which without having to get back to the actual provider website.

The third step in defining a dataset provider is to overload, for each
`LayerProvider` and `LayerDataset`, the `_rasterpath` method; this method is
responsible for creating a path where the layers themselves (and when
applicable, the relevant zip/tar files) will be stored, relative to the root of
the data storage folder. Note that by default, there is a method that takes two
arguments (the provider and the dataset) and returns the joined paths of each;
you can overload this two-argument methods if you want the same dataset to be
stored in differently named folders for each provider. Note also that the
two-method argument does not exists for scenarios (for the moment). Note finally
that these functions do not allow *users* to determine where the layers will be
stored, as they are meant to access them using overloads of `SimpleSDMPredictor`
anyways.