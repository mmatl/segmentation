"""
File for picking unique, easily-distinguishable colors.
Author: Matthew Matl
"""

def indexed_color(i):
    """Return a color based on an index. Identical indices will always return
    the same color, and the colors should contrast cleanly.

    Parameters
    ----------
    i : int
        An index into the color array.

    Returns
    -------
    :obj:`tuple` of float
        An rgb color tuple.
    """

    r = [0.95, 0.13, 0.95, 0.53, 0.95, 0.63, 0.75, 0.76, 0.52, 0.0,
            0.9, 0.0, 0.98, 0.38, 0.96, 0.7, 0.86, 0.53, 0.55, 0.4, 0.89, 0.17]
    g = [0.95, 0.13, 0.76, 0.34, 0.52, 0.79, 0.0, 0.7, 0.52, 0.53,
            0.56, 0.4, 0.58, 0.31, 0.65, 0.27, 0.83, 0.18, 0.71, 0.27, 0.35, 0.24]
    b = [0.96, 0.13, 0.0, 0.57, 0.0, 0.95, 0.2, 0.5, 0.51, 0.34, 0.67,
            0.65, 0.47, 0.59, 0.0, 0.42, 0.0, 0.09, 0.0, 0.13, 0.13, 0.15]
    ci = i % len(r)
    return (r[ci], g[ci], b[ci])
