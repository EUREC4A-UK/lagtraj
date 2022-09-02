import datetime

import isodate
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
            yield (prefix, item)
        elif isinstance(item, float) or isinstance(item, int):
            yield (prefix, str(item))
        elif isinstance(item, datetime.timedelta):
            yield (prefix, isodate.duration_isoformat(item))
        elif isinstance(item, datetime.datetime):
            yield (prefix, item.isoformat())
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
            # can use isinstance because namedtuple is a kind of tuple and we
            # want to save the tuple keys in that case
            elif isinstance(item, list) or type(item) == tuple:
                sub_items = zip(range(len(item)), item)
            else:
                # collections.named_tuple has a `_asdict` method to turn it into a dictionary
                if hasattr(item, "_asdict"):
                    sub_items = item._asdict().items()
                else:
                    sub_items = filter(
                        lambda item: not item[0].startswith("_"), vars(item).items()
                    )

            for (k, v) in sub_items:
                # skip None values
                if v is None:
                    continue

                if prefix != "":
                    label = f"{prefix}_{k}"
                else:
                    label = k

                for serialized_item in _serialize_item(v, prefix=label):
                    yield serialized_item

    attrs = {}
    items = list(args) + [
        kwargs,
    ]
    for item in items:
        for (label, value) in _serialize_item(item=item):
            attrs[label] = value
    return attrs
