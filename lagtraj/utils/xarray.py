import xarray as xr
from ..domain.sources.era5.load import ERA5DataSet


def create_attributes_dictionary(*args, **kwargs):
    """
    Serialises a (possibly nested) dictionary structure into a one-level
    dictionary by joining the keys at each level with underscore ("_"). For
    xarray Datasets and DataArrays the attributes are used.
    """

    def _serialize_item(item, prefix=""):
        if type(item) == str:
            yield (label, item)
        elif type(item) in [float, int]:
            yield (label, str(item))
        else:
            sub_items = []
            if type(item) == dict:
                sub_items = item.items()
            elif (
                isinstance(item, xr.Dataset)
                or isinstance(item, xr.DataArray)
                or isinstance(item, ERA5DataSet)
            ):
                sub_items = item.attrs.items()
            else:
                sub_items = filter(
                    lambda item: not item[0].startswith("_"), vars(item).items()
                )

            for (k, v) in items:
                if prefix != "":
                    label = f"{prefix}_{k}"
                else:
                    label = k
                _serialize_item(v, prefix=label)

    attrs = {}
    items = list(args) + [
        kwargs,
    ]
    for item in items:
        for (label, value) in _serialize_item(item=item):
            attrs[label] = value
    return attrs
