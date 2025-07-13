#!/usr/bin/env python3
# This script creates a set of test points for WildF functions.

from pathlib import Path
import numpy as np
import mpmath
from math import log10
import os

os.makedirs('test/test_points', exist_ok=True)
os.makedirs('test/test_specfun', exist_ok=True)

mpmath.mp.dps = 16
test_points_dir = Path(__file__).parent / 'test_points/'

def set_test_points_exponential_integral_ei():
    # Ei(x) goes to infinity as x increases. Can't test it for x > 709.0.
    xc = np.logspace(-100, -1, 100, dtype=np.float64)
    xf = np.linspace(1.0, 709.0, 1000, dtype=np.float64)
    xv = np.concatenate((xc, xf))
    yv = np.array([mpmath.ei(x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'exponential_integral_ei.csv'
    np.savetxt(filepath, zv, fmt='% .15e')

def set_test_points_exponential_integral_e1x():
    # E1(x) goes to 0 as x goes to infinity. No need to test for x > 10³.
    xc = np.logspace(-100, -1, 100, dtype=np.float64)
    xf = np.linspace(1.0, 739.0, 1000, dtype=np.float64)
    xv = np.concatenate((xc, xf))
    yv = np.array([mpmath.e1(x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'exponential_integral_e1x.csv'
    np.savetxt(filepath, zv, fmt='% .15e')

def set_test_points_exponential_integral_e1z():
    # {z ∈ ℂ | Re(z) ≥ 0, z ≠ 0}
    x1c = np.logspace(-100, -1, 25, dtype=np.float64)
    x1f = np.linspace(1.0, 900.0, 25, dtype=np.float64)
    x1w = np.logspace(3, 16, 14, dtype=np.float64)
    x1 = np.concatenate(([0.0], x1c, x1f, x1w), dtype=np.float64)
    y1c = np.logspace(-100, -1, 25, dtype=np.float64)
    y1f = np.linspace(1.0, 900.0, 25, dtype=np.float64)
    y1w = np.logspace(3, 16, 14, dtype=np.float64)
    y1p = np.concatenate(([0.0], y1c, y1f, y1w), dtype=np.float64)
    y1n = -y1p[::-1].copy()
    y1 = np.concatenate((y1n, [0.0], y1p), dtype=np.float64)
    U1, V1 = np.meshgrid(x1, y1)
    W1 = U1 + 1j * V1
    z1 = W1.ravel()
    z1_mask = ~np.logical_and(z1.real == 0, z1.imag == 0)
    z1 = z1[z1_mask]

    # {z ∈ ℂ | Re(z) < 0, |Im(z)| ≥ 1e-6, Im(z) = 0}
    x2c = np.logspace(-100, -1, 25, dtype=np.float64)
    x2f = np.linspace(1.0, 700.0, 25, dtype=np.float64)
    x2 = -np.concatenate((x2c, x2f))[::-1]
    y2p = np.logspace(-6, 8, 15, dtype=np.float64)
    y2n = -y2p[::-1].copy()
    y2 = np.concatenate((y2n, [0.0], y2p), dtype=np.float64)
    U2, V2 = np.meshgrid(x2, y2)
    W2 = U2 + 1j * V2
    z2 = W2.ravel()

    zv = np.concatenate((z1, z2))
    wv = np.array([mpmath.e1(z) for z in zv], dtype=np.complex128)
    uv = np.column_stack((zv, wv))
    filepath = test_points_dir / 'exponential_integral_e1z.csv'
    np.savetxt(filepath, uv, fmt='% .15e %+.15e   % .15e %+.15e')

def set_test_points_bessel_j0x():
    # J0(x) goes to 0 as x goes to infinity. No need to test for x > 1e32.
    xp = np.logspace(-8, 32, 100, dtype=np.float64)
    xv = np.concatenate(([0.0], xp), dtype=np.float64)
    yv = np.array([mpmath.besselj(0, x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'bessel_j0x.csv'
    np.savetxt(filepath, zv, fmt='% .15e')

def set_test_points_bessel_j1x():
    # J1(x) goes to 0 as x goes to infinity. No need to test for x > 1e32.
    xp = np.logspace(-8, 32, 100, dtype=np.float64)
    xv = np.concatenate(([0.0], xp), dtype=np.float64)
    yv = np.array([mpmath.besselj(1, x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'bessel_j1x.csv'
    np.savetxt(filepath, zv, fmt='% .15e')

def set_test_points_bessel_y0x():
    # Y0(x) goes to 0 as x goes to infinity. No need to test for x > 1e32.
    xv = np.logspace(-17, 32, 100, dtype=np.float64)
    yv = np.array([mpmath.bessely(0, x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'bessel_y0x.csv'
    np.savetxt(filepath, zv, fmt='% .15e')

def set_test_points_bessel_y1x():
    # Y1(x) goes to 0 as x goes to infinity. No need to test for x > 1e32.
    xv = np.logspace(-17, 32, 100, dtype=np.float64)
    yv = np.array([mpmath.bessely(1, x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'bessel_y1x.csv'
    np.savetxt(filepath, zv, fmt='% .15e')

def set_test_points_hypergeometric_hyp2f1():
    # In the test for hyp2f1(a,b,c,z), all parameters are of type complex with zero
    # imaginary part. Parameters a anc c are fixed. Parameter b ranges through
    # multiples of -n/2, where n is a non-negative integer. Parameter z ranges through
    # values larger or equal to 1.0.
    a = 0.5 + 0.0j
    c = 1.5 + 0.0j
    
    n = np.arange(0, 42, dtype=np.complex128)
    b = -0.5 * n[::-1]

    z = np.linspace(1.0, 5.0, 100, dtype=np.complex128)

    B, Z = np.meshgrid(b, z)

    bv = B.ravel()
    zv = Z.ravel()
    hv = np.empty_like(bv)
    for i, (bi, zi) in enumerate(zip(bv, zv)):
        hv[i] = mpmath.hyp2f1(a, bi, c, zi)

    xv = bv + 1j * zv  # A subterfuge to allow fortran test routines to read this data.
    dv = np.column_stack((xv, hv))
    filepath = test_points_dir / 'hypergeometric_hyp2f1.csv'
    np.savetxt(filepath, dv, fmt='% .15e %+.15e   % .15e %+.15e')


if __name__ == '__main__':
    set_test_points_exponential_integral_ei()
    set_test_points_exponential_integral_e1x()
    set_test_points_exponential_integral_e1z()
    set_test_points_bessel_j0x()
    set_test_points_bessel_j1x()
    set_test_points_bessel_y0x()
    set_test_points_bessel_y1x()
    set_test_points_hypergeometric_hyp2f1()
