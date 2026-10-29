"""NumPy fallback for the optional compiled KMC kernels."""

import numpy as np


def distances(dx, dy, dz, num):
    dx = np.asarray(dx, dtype=float)
    dy = np.asarray(dy, dtype=float)
    dz = np.asarray(dz, dtype=float)
    return np.sqrt(dx[:num] ** 2 + dy[:num] ** 2 + dz[:num] ** 2)


def forster(Rf, mats, num, alpha_mu, r, emi_rate):
    r = np.asarray(r, dtype=float)[:num]
    mats = np.asarray(mats, dtype=int)[:num]
    rates = np.zeros(num, dtype=float)
    mask = r != 0
    ratio = np.asarray(Rf, dtype=float)[mats[mask]] / (alpha_mu + r[mask])
    rates[mask] = ratio ** 6 * emi_rate
    return rates


def forster_anni(
    Rf,
    mats,
    num,
    alpha_mu,
    r,
    emi_rate,
    replace_pos,
    replace_raios,
    mum,
):
    rates = forster(Rf, mats, num, alpha_mu, r, emi_rate)
    positions = np.asarray(replace_pos, dtype=int)[:mum]
    radii = np.asarray(replace_raios, dtype=float)[:mum]
    distances_ = np.asarray(r, dtype=float)[positions]
    rates[positions] = (radii / (alpha_mu + distances_)) ** 6 * emi_rate
    return rates


def jump(jump_rate, num, random_number):
    rates = np.asarray(jump_rate, dtype=float)[:num]
    total = float(rates.sum())
    if num == 0 or total <= 0:
        return 0.0, -1
    index = int(
        np.searchsorted(np.cumsum(rates), random_number * total, side="left")
    )
    return total, min(index, num - 1)


def dexter(Rd, invL, emi_rate, mats, num, r):
    r = np.asarray(r, dtype=float)[:num]
    mats = np.asarray(mats, dtype=int)[:num]
    radii = np.asarray(Rd, dtype=float)[mats]
    rates = np.zeros(num, dtype=float)
    mask = r != 0
    rates[mask] = (
        (radii[mask] * invL) ** 2
        * np.exp(2 * invL * (radii[mask] - r[mask]))
        * emi_rate
    )
    return rates
