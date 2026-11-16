"""
Structural integrity: habitat quality of the pixel and of its surrounding landscape.

Structural integrity follows Hill et al. (2026), who apply the landscape quality of Beyer et
al. (2020) to human modification:

    S_i = sqrt(w_i) * sum_j K_ij sqrt(w'_j) / sum_j K_ij,   K_ij = exp(-d_ij / 5 km)

with the kernel cut at 26.5 km (``LANDSCAPE_DECAY_M``, ``LANDSCAPE_RADIUS_M``). Habitat quality
``w = habitat_quality(HMI)`` is a logistic curve with midpoint HMI 0.4 (``QUALITY_MIDPOINT``,
calibrated against WorldCover natural land cover). The site term uses ``w``; the landscape term
uses ``w'``, the grey opening of ``w`` by a 300 m ground disk (``LANDSCAPE_OPENING_M``), so
isolated patches narrower than about 600 m add no landscape quality. Uniform land of HMI h
scores ``habitat_quality(h)``. There is no habitat threshold and no class step.

Kernels are ellipses in pixel space that are circles on the ground: their east-west radius
grows with 1/cos(latitude). Pixels without HMI (sea, nodata) do not enter the kernel mean and
stay NaN in the output.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.ndimage import grey_opening
from scipy.signal import fftconvolve

from .settings import (
    LANDSCAPE_DECAY_M,
    LANDSCAPE_OPENING_M,
    LANDSCAPE_RADIUS_M,
    METERS_PER_DEGREE,
    QUALITY_MIDPOINT,
    QUALITY_STEEPNESS,
)

STRIP_DEG = 0.5
"""Rows within this latitude span share one kernel (radius error < 1 % below 75 deg)."""


def habitat_quality(
    hmi: np.ndarray | float,
    midpoint: float = QUALITY_MIDPOINT,
    steepness: float = QUALITY_STEEPNESS,
) -> np.ndarray:
    """Habitat quality of HMI: 1 at HMI 0, 0 at HMI 1.

    The logistic ``1 / (1 + exp(steepness * (HMI - midpoint)))``, rescaled to span [0, 1] on
    HMI [0, 1]. Hill et al. (2026) use the linear ``1 - HMI``.

    Example:
        >>> habitat_quality(np.array([0.1, 0.3, 0.5, 0.7])).round(2)
        array([0.97, 0.74, 0.27, 0.05])
    """
    hmi = np.asarray(hmi, dtype="float64")

    def logistic(h: np.ndarray | float) -> np.ndarray:
        return 1 / (1 + np.exp(steepness * (np.asarray(h) - midpoint)))

    return np.asarray((logistic(hmi) - logistic(1.0)) / (logistic(0.0) - logistic(1.0)))


def ground_disk(radius_m: float, lat: float, pixel_deg: float) -> np.ndarray:
    """Boolean footprint of a ground circle on a lat/lon grid at latitude ``lat``.

    Like Earth Engine's ``Kernel.circle(units="meters")``, a pixel belongs to the kernel
    when its centre is within the radius, so 300 m on a 300 m grid is a 3x3 cross.
    """
    ry = radius_m / (pixel_deg * METERS_PER_DEGREE) + 1e-9
    rx = ry / max(math.cos(math.radians(lat)), 0.01)
    y, x = np.ogrid[-math.floor(ry) : math.floor(ry) + 1, -math.floor(rx) : math.floor(rx) + 1]
    return (y / ry) ** 2 + (x / rx) ** 2 <= 1


def exponential_kernel(radius_m: float, decay_m: float, lat: float, pixel_deg: float) -> np.ndarray:
    """``exp(-d / decay_m)`` on the pixels of ``ground_disk(radius_m, lat, pixel_deg)``.

    Ground distances use the local scale, ``cos(lat)`` metres east-west per north-south metre.
    """
    disk = ground_disk(radius_m, lat, pixel_deg)
    ny, nx = disk.shape[0] // 2, disk.shape[1] // 2
    y, x = np.ogrid[-ny : ny + 1, -nx : nx + 1]
    scale = max(math.cos(math.radians(lat)), 0.01)
    d = pixel_deg * METERS_PER_DEGREE * np.hypot(y, x * scale)
    return np.where(disk, np.exp(-d / decay_m), 0).astype("float32")


def structural_halo_m() -> float:
    """Ground distance a window must extend beyond the pixels whose integrity is used."""
    return LANDSCAPE_RADIUS_M + LANDSCAPE_OPENING_M


def structural_integrity(hmi: np.ndarray, top_lat: float, pixel_deg: float) -> np.ndarray:
    """Structural integrity of every pixel of an HMI window on an EPSG:4326 grid.

    Pixels near the window edge only see part of their neighbourhood; read the window with
    a halo of ``structural_halo_m()`` and crop.

    Args:
        hmi: Human Modification Index (0-1), NaN = no data, square pixels of ``pixel_deg``.
        top_lat: Latitude of the window's top edge.
        pixel_deg: Pixel size in degrees.

    Returns:
        float32 structural integrity in [0, 1], NaN where HMI is NaN.
    """
    return landscape_quality(hmi, top_lat, pixel_deg)


def landscape_quality(
    hmi: np.ndarray,
    top_lat: float,
    pixel_deg: float,
    decay_m: float = LANDSCAPE_DECAY_M,
    radius_m: float = LANDSCAPE_RADIUS_M,
    opening_m: float = LANDSCAPE_OPENING_M,
    midpoint: float = QUALITY_MIDPOINT,
) -> np.ndarray:
    """Structural integrity, ``sqrt(w_i)`` times the kernel mean of ``sqrt(opened w)``.

    Args:
        hmi: Human Modification Index (0-1), NaN = no data, square pixels of ``pixel_deg``.
        top_lat: Latitude of the window's top edge.
        pixel_deg: Pixel size in degrees.
        decay_m: Exponential kernel scale.
        radius_m: Kernel truncation radius.
        opening_m: Radius of the grey opening of the landscape term (0 = none).
        midpoint: Midpoint of ``habitat_quality``.

    Returns:
        float32 structural integrity in [0, 1], NaN where HMI is NaN.
    """
    valid = np.isfinite(hmi)
    # No-data carries w = 1 so that sea does not open away coastal land; it never enters
    # the kernel mean.
    quality = habitat_quality(np.clip(np.nan_to_num(hmi, nan=0.0), 0, 1), midpoint)
    w = np.where(valid, quality, 1).astype("float32")
    out = np.full(hmi.shape, np.nan, dtype="float32")
    rows = hmi.shape[0]
    strip = max(1, round(STRIP_DEG / pixel_deg))
    for r0 in range(0, rows, strip):
        r1 = min(r0 + strip, rows)
        lat = top_lat - (r0 + r1) / 2 * pixel_deg
        kernel = exponential_kernel(radius_m, decay_m, lat, pixel_deg)
        footprint = ground_disk(opening_m, lat, pixel_deg) if opening_m > 0 else None
        # Rows the mean of r0:r1 needs: the kernel radius, plus erosion and dilation below it.
        pad = kernel.shape[0] // 2 + (0 if footprint is None else 2 * (footprint.shape[0] // 2))
        a0, a1 = max(r0 - pad, 0), min(r1 + pad, rows)
        land = w[a0:a1]
        if footprint is not None:
            land = grey_opening(land, footprint=footprint, mode="nearest")
        v = valid[a0:a1].astype("float32")
        total = fftconvolve(np.sqrt(land) * v, kernel, mode="same")
        count = fftconvolve(v, kernel, mode="same")
        with np.errstate(divide="ignore", invalid="ignore"):
            mean = total / count
        out[r0:r1] = np.clip(np.sqrt(w[r0:r1]) * mean[r0 - a0 : r1 - a0], 0, 1)
    out[~valid] = np.nan
    return out
