#!/usr/bin/env python3
# This script creates a markdown file with the test results.

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import os

os.makedirs('test/test_plots', exist_ok=True)

test_specfun_dir = Path(__file__).parent / 'test_specfun/'
test_plots_dir = Path(__file__).parent / 'test_plots/'

raw_url = 'https://raw.githubusercontent.com/rodpcastro/colspecf/refs/heads/gauss/test/test_plots/'

eps = np.finfo(np.float64).eps


def float_conv(s):
    new_s = s
    if 'E' not in s:
        s1 = s[1:]
        esign = '+' if '+' in s1 else '-'
        s_split = s1.split(esign)
        new_s = s[0] + s_split[0] + 'E' + esign + s_split[1]

    return np.float64(new_s)


def compute_error(ym, yc):
    if ym.dtype == np.float64:
        abs_rel_err = np.empty((2, *ym.shape))
        abs_rel_err[0] = np.abs(yc - ym)
        abs_rel_err[1] = abs_rel_err[0] / (np.abs(ym) + eps)
        err = np.min(abs_rel_err, axis=0)
    elif ym.dtype == np.complex128:
        err_re_im = np.empty((2, *ym.shape))
        err_re_im[0] = compute_error(ym.real, yc.real)
        err_re_im[1] = compute_error(ym.imag, yc.imag)
        err = np.max(err_re_im, axis=0)

    return err


def write_results(f, results_title, get_error, fname, latex):
    err = get_error(fname, latex)
    f.write(r'## ' + results_title + '\n\n')
    f.write(f'| Maximum Error | Average Error  |\n')
    f.write(f'| ------------- | -------------- |\n')
    f.write(f'|{err.max():.2e}|{err.mean():.2e}|\n\n')
    f.write(f'![{fname}]({raw_url + fname + '.svg'})\n\n')


def get_error_fx(filename, latex=''):
    filepath = test_specfun_dir / (filename + '.csv')
    
    columns = list(range(3))
    converters = dict([(i, float_conv) for i in columns])
    
    data = np.loadtxt(filepath, usecols=columns, converters=converters)
    
    x = data[:, 0]
    ym = data[:, 1]
    yc = data[:, 2]
    
    err = compute_error(ym, yc)
    
    fig, ax = plt.subplots(figsize=(5, 3))
    
    ax.set_title(rf'Error - {latex}')
    ax.set_xlabel('x')
    ax.plot(x, err, 'k')

    imgpath = test_plots_dir / (filename + '.svg')
    plt.savefig(imgpath, bbox_inches='tight')
    
    return err


def get_error_fz(filename, latex=''):
    filepath = test_specfun_dir / (filename + '.csv')
    
    columns = list(range(6))
    converters = dict([(i, float_conv) for i in columns])
    
    data = np.loadtxt(filepath, usecols=columns, converters=converters)
    
    n = len(np.where(data[:, 0] == data[0, 0])[0])
    m = len(np.where(data[:, 1] == data[0, 1])[0])
    
    zre = data[:, 0].reshape((n, m))
    zim = data[:, 1].reshape((n, m))
    ym = data[:, 2] + 1j * data[:, 3]
    ym = ym.reshape((n, m))
    yc = data[:, 4] + 1j * data[:, 5]
    yc = yc.reshape((n, m))
    
    err = compute_error(ym, yc)
    
    fig, ax = plt.subplots(figsize=(5, 3))
    
    ax.set_title(rf'Error - {latex}')
    ax.set_xlabel(r'$\Re(z)$')
    ax.set_ylabel(r'$\Im(z)$', rotation='horizontal', labelpad=10.0)
    ec = ax.contourf(zre, zim, err)
    
    fig.colorbar(ec, format='%.2e')
    
    imgpath = test_plots_dir / (filename + '.svg')
    plt.savefig(imgpath, bbox_inches='tight')
    
    return err


def get_error_hypergeometric_hyp2f1(filename, latex):
    filepath = test_specfun_dir / (filename + '.csv')

    columns = list(range(6))
    converters = dict([(i, float_conv) for i in columns])

    data = np.loadtxt(filepath, usecols=columns, converters=converters)

    n = len(np.where(data[:, 0] == data[0, 0])[0])
    m = len(np.where(data[:, 1] == data[0, 1])[0])

    b = data[:, 0].reshape((n, m))
    z = data[:, 1].reshape((n, m))
    hm = data[:, 2] + 1j * data[:, 3]
    hm = hm.reshape((n, m))
    hc = data[:, 4] + 1j * data[:, 5]
    hc = hc.reshape((n, m))

    err = compute_error(hm, hc)

    fig, ax = plt.subplots(figsize=(5, 3))

    ax.set_title(rf'Error - {latex}')
    ax.set_xlabel('b')
    ax.set_ylabel('z', rotation='horizontal', labelpad=10.0)
    ec = ax.contourf(b, z, err)
    
    fig.colorbar(ec, format='%.2e')
    
    imgpath = test_plots_dir / (filename + '.svg')
    plt.savefig(imgpath, bbox_inches='tight')
    
    return err


if __name__ == '__main__':
    filepath = Path(__file__).parent / 'test_results.md'
    with open(filepath, 'w') as f:
        f.write('# Test results\n\n')
        latex = r'$\mathrm{Ei}(x)$'
        title = rf'Exponential integral {latex}'
        fname = 'exponential_integral_ei'
        write_results(f, title, get_error_fx, fname, latex)
    
        latex = r'$\mathrm{E}_1(x)$'
        title = rf'Exponential integral {latex}'
        fname = 'exponential_integral_e1x'
        write_results(f, title, get_error_fx, fname, latex)
    
        latex = r'$\mathrm{E}_1(z)$'
        title = rf'Exponential integral {latex}'
        fname = 'exponential_integral_e1z'
        write_results(f, title, get_error_fz, fname, latex)
        
        latex = '$J_0(x)$'
        title = rf'Bessel function {latex}'
        fname = 'bessel_j0x'
        write_results(f, title, get_error_fx, fname, latex)
        
        latex = '$J_1(x)$'
        title = rf'Bessel function {latex}'
        fname = 'bessel_j1x'
        write_results(f, title, get_error_fx, fname, latex)
        
        latex = '$Y_0(x)$'
        title = rf'Bessel function {latex}'
        fname = 'bessel_y0x'
        write_results(f, title, get_error_fx, fname, latex)
        
        latex = '$Y_1(x)$'
        title = rf'Bessel function {latex}'
        fname = 'bessel_y1x'
        write_results(f, title, get_error_fx, fname, latex)
    
        latex = '${}_2F_1(a, b; c; z)$'
        title = rf'Gauss hypergeometric function {latex}'
        fname = 'hypergeometric_hyp2f1'
        write_results(f, title, get_error_hypergeometric_hyp2f1, fname, latex)

