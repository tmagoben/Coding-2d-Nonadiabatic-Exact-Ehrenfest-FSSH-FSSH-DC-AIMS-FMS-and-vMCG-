"""Uniform periodic grids and their Fourier conjugate momenta."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class PeriodicGrid1D:
    """A uniform periodic coordinate grid compatible with the FFT."""

    x: FloatArray
    k: FloatArray
    dx: float
    length: float

    @classmethod
    def create(cls, n: int, x_min: float, x_max: float) -> "PeriodicGrid1D":
        if n < 4:
            raise ValueError("A useful periodic grid requires at least 4 points.")
        if x_max <= x_min:
            raise ValueError("x_max must be greater than x_min.")

        length = float(x_max - x_min)
        dx = length / n
        x = x_min + dx * np.arange(n, dtype=float)
        k = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
        return cls(x=x, k=k, dx=dx, length=length)


@dataclass(frozen=True)
class PeriodicGrid2D:
    """A tensor-product periodic grid compatible with two-dimensional FFTs."""

    x: FloatArray
    y: FloatArray
    kx: FloatArray
    ky: FloatArray
    dx: float
    dy: float
    length_x: float
    length_y: float

    @classmethod
    def create(
        cls,
        nx: int,
        ny: int,
        x_min: float,
        x_max: float,
        y_min: float,
        y_max: float,
    ) -> "PeriodicGrid2D":
        if nx < 4 or ny < 4:
            raise ValueError("Each grid direction requires at least 4 points.")
        if x_max <= x_min or y_max <= y_min:
            raise ValueError("Grid maxima must be greater than grid minima.")

        length_x = float(x_max - x_min)
        length_y = float(y_max - y_min)
        dx = length_x / nx
        dy = length_y / ny

        x = x_min + dx * np.arange(nx, dtype=float)
        y = y_min + dy * np.arange(ny, dtype=float)
        kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
        ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy)

        return cls(
            x=x,
            y=y,
            kx=kx,
            ky=ky,
            dx=dx,
            dy=dy,
            length_x=length_x,
            length_y=length_y,
        )
