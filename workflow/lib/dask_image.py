from __future__ import annotations

import os
from glob import glob

try:
    import tifffile
except (AttributeError, ImportError):
    pass

from dask.array.core import Array
from dask.base import tokenize


def add_leading_dimension(x):
    return x[None, ...]

def imread(fn, page):
    return tifffile.imread(fn,key=page)


def imread_pages(filename, preprocess=None):

    
    tif = tifffile.TiffFile(filename)
    pages = [i for i in range(len(tif.pages))]
    name = "imread-%s" % tokenize(filename, os.path.getmtime(filename))

    sample = tif.pages[0].asarray()

    if preprocess:
        sample = preprocess(sample)

    keys = [(name, i) + (0,) * len(sample.shape) for i in pages]
    if preprocess:
        values = [
            (add_leading_dimension, (preprocess, (imread, filename,i))) for i in pages
        ]
    else:
        values = [(add_leading_dimension, (imread, filename,i)) for i in pages]
    dsk = dict(zip(keys, values))

    chunks = ((1,) * len(pages),) + tuple((d,) for d in sample.shape)

    return Array(dsk, name, chunks, sample.dtype)
