#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Util code - some methods from https://github.com/cerinunn/pdart
This project only needs these three methods. By copying them,
we don't need any dependencies.
Author: Ceri Nunn, JPL
"""
import numpy as np
import numpy.ma as ma
import pandas as pd

def remove_negative_ones(stream,channels=['MH1','MH2','MHZ','SHZ','ATT']):
    """Snippet to remove the -1 values in the data traces.

    The SHZ traces have missing data samples 3-4 times every 32 samples.
    Providing the seed data with these missing data would mean using very
    large files. Instead, we provide the data with -1 replacing the gaps.
    To change the files to include the gaps, use this simple method to
    replace the -1 values.
    """

    for tr in stream:
        if tr.stats.channel in channels:
            if tr.stats.channel in ('MH1','MH2','MHZ','SHZ'):
                tr.data = np.ma.masked_where(tr.data==-1, tr.data)
            elif tr.stats.channel in ('ATT'):
                tr.data = np.ma.masked_values(tr.data,-1.0)

def remove_negative_ones_trace(trace):
    """Snippet to remove the -1 values in the data traces.

    The SHZ traces have missing data samples 3-4 times every 32 samples.
    Providing the seed data with these missing data would mean using very
    large files. Instead, we provide the data with -1 replacing the gaps.
    To change the files to include the gaps, use this simple method to
    replace the -1 values.
    """

    if trace.stats.channel in ('MH1','MH2','MHZ','SHZ'):
        trace.data = np.ma.masked_where(trace.data==-1, trace.data)
    elif trace.stats.channel in ('ATT'):
        trace.data = np.ma.masked_values(trace.data,-1.0)

def linear_interpolation(trace,interpolation_limit=1):
    """Snippet to interpolate missing data.

    The SHZ traces have missing data samples 3-4 times every 32 samples.
    Providing the seed data with these missing data would mean using very
    large files. Instead, we provide the data with -1 replacing the gaps.
    To change the files to interpolate across the gaps, use this simple method to
    replace the -1 values. The trace is modified, and a mask is applied at
    the end if necessary.

    :type stream: :class:`~obspy.core.Trace`
    :param trace: A data trace
    :type interpolation_limit: int
    :param interpolation_limit: Limit for interpolation. Defaults to 1. For
      more information read the options for the `~pandas.Series.interpolate`
      method.

    :return: original_mask :class:`~numpy.ndarray` or class:`~numpy.bool_`
       Returns the original mask, before any interpolation is made.

    """

    trace.data = np.ma.masked_where(trace.data == -1, trace.data)
    original_mask = np.ma.getmask(trace.data)
    data_series = pd.Series(trace.data)
    # data_series.replace(-1.0, pd.NA, inplace=True)
    data_series.interpolate(method='linear', axis=0, limit=interpolation_limit, inplace=True, limit_direction=None, limit_area='inside', downcast=None)
    data_series.fillna(-1.0, inplace=True)
    trace.data=data_series.to_numpy(dtype=int)
    trace.data = np.ma.masked_where(trace.data == -1, trace.data)
    return original_mask
