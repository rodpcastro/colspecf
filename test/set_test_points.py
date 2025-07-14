#!/usr/bin/env python3
# This script creates a set of test points for ColSpecF functions.

from pathlib import Path
import numpy as np
import mpmath
import os

os.makedirs('test/test_points', exist_ok=True)
os.makedirs('test/test_specfun', exist_ok=True)

mpmath.mp.dps = 16
test_points_dir = Path(__file__).parent / 'test_points/'


def set_test_points_exponential_integral_ei():
    # Ei(x) goes to infinity as x increases. Can't test it for x > 709.0.
    xv = np.linspace(1.0e-8, 100.0, 100)
    yv = np.array([mpmath.ei(x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'exponential_integral_ei.csv'
    np.savetxt(filepath, zv, fmt='% .15e')


def set_test_points_exponential_integral_e1x():
    # E1(x) goes to 0 as x goes to infinity. No need to test for x > 10³.
    xv = np.linspace(1.0e-8, 100.0, 100)
    yv = np.array([mpmath.e1(x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'exponential_integral_e1x.csv'
    np.savetxt(filepath, zv, fmt='% .15e')


def set_test_points_exponential_integral_e1z():
    x = np.linspace(-100.0, 100.0, 100)
    y = np.linspace(-100.0, 100.0, 100)
    X, Y = np.meshgrid(x, y) 
    Z = X + 1j * Y
    zv = Z.ravel()
    wv = np.array([mpmath.e1(z) for z in zv], dtype=np.complex128)
    uv = np.column_stack((zv, wv))
    filepath = test_points_dir / 'exponential_integral_e1z.csv'
    np.savetxt(filepath, uv, fmt='% .15e %+.15e   % .15e %+.15e')


def set_test_points_bessel_j0x():
    # J0(x) goes to 0 as x goes to infinity. No need to test for x > 1e32.
    xv = np.linspace(0.0, 100.0, 100)
    yv = np.array([mpmath.besselj(0, x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'bessel_j0x.csv'
    np.savetxt(filepath, zv, fmt='% .15e')


def set_test_points_bessel_j1x():
    # J1(x) goes to 0 as x goes to infinity. No need to test for x > 1e32.
    xv = np.linspace(0.0, 100.0, 100)
    yv = np.array([mpmath.besselj(1, x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'bessel_j1x.csv'
    np.savetxt(filepath, zv, fmt='% .15e')

def set_test_points_bessel_y0x():
    # Y0(x) goes to 0 as x goes to infinity. No need to test for x > 1e32.
    xv = np.linspace(1.0e-8, 100.0, 100)
    yv = np.array([mpmath.bessely(0, x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'bessel_y0x.csv'
    np.savetxt(filepath, zv, fmt='% .15e')

def set_test_points_bessel_y1x():
    # Y1(x) goes to 0 as x goes to infinity. No need to test for x > 1e32.
    xv = np.linspace(1.0e-8, 100.0, 100)
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
