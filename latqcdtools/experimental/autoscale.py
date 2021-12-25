import numpy as np


def auto_range(data, accept=0.3):
    data_new = []
    for value in data:
        if not (np.isnan(value) or np.isinf(value)):
            data_new.append(value)
    data = data_new
    mean = np.median(data)
    mn = mean
    mx = mean
    for i in data:
        if mean > 0:
            if (1 - accept)*mean < i < (1 + accept)*mean:
                mn = min(mn, i)
                mx = max(mx, i)
        else:
            if (1 - accept)*mean > i > (1 + accept)*mean:
                mn = min(mn, i)
                mx = max(mx, i)
    return mn, mx


def auto_range_multi(data, accept=0.3):
    mn, mx = auto_range(data[0], accept)
    for i in data[1:]:
        mn_new, mx_new = auto_range(i, accept)
        mn = min(mn, mn_new)
        mx = max(mx, mx_new)
    return mn, mx


def auto_range_err(data, data_err, accept=0.3):
    data = np.array(data)
    data_err = np.array(data_err)
    tmp = [data - data_err, data + data_err]
    return auto_range_multi(tmp, accept)
