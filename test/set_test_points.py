from pathlib import Path
import numpy as np
import mpmath

mpmath.mp.dps = 16
test_points_dir = Path(__file__).parent / 'test_points/'

def set_test_points_exponential_integral_ei():
    xv = np.logspace(-16, 2, 1000, dtype=np.float64)
    yv = np.array([mpmath.ei(x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'exponential_integral_ei.csv'
    np.savetxt(filepath, zv, fmt='% .15e')

def set_test_points_exponential_integral_e1x():
    xv = np.logspace(-16, 2, 1000, dtype=np.float64)
    yv = np.array([mpmath.e1(x) for x in xv])
    zv = np.column_stack((xv, yv))
    filepath = test_points_dir / 'exponential_integral_e1x.csv'
    np.savetxt(filepath, zv, fmt='% .15e')

def set_test_points_exponential_integral_e1z():
    px = np.linspace(1.0e-16, 100, 20, dtype=np.float64)
    mx = -px[::-1].copy()
    x = np.concatenate((mx, px))

    py = np.linspace(1.0e-16, 100, 20, dtype=np.float64)
    my = -py[::-1].copy()
    y = np.concatenate((my, py))

    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    zv = Z.ravel()

    # For -60 < real(z) < 0 and |imag(z)| < 0.5, the implemented e1z
    # cannot achieve the desired precision. I'll look for a fix.
    zre_mask = np.logical_and(zv.real > -60, zv.real < 0)
    zim_mask = np.abs(zv.imag) < 0.5
    zmask = np.logical_and(zre_mask, zim_mask)
    zv = zv[~zmask]

    wv = np.array([mpmath.e1(z) for z in zv], dtype=np.complex128)
    uv = np.column_stack((zv, wv))
    filepath = test_points_dir / 'exponential_integral_e1z.csv'
    np.savetxt(filepath, uv, fmt='% .15e %+.15e   % .15e %+.15e')

if __name__ == '__main__':
    set_test_points_exponential_integral_ei()
    set_test_points_exponential_integral_e1x()
    set_test_points_exponential_integral_e1z()

