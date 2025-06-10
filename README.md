<p align="center">
    <img src="https://raw.githubusercontent.com/rodpcastro/colspecf/refs/heads/main/media/colspecf_logo.svg", alt="ColSpecF Logo">
</p>

[![Fortran][Fortran Badge]][Fortran Website]
[![Docs][Docs Badge]][Docs Website]
[![Test][Test Badge]][Test Workflow]
[![Coverage][Coverage Badge]][Coverage Website]
[![License][License Badge]][License File]

# ColSpecF
ColSpecF (Collected Special Functions) is a [Fortran][Fortran Website] library for evaluating mathematical [Special Functions], built around adaptations of [Collected Algorithms][calgo] from [ACM] to modern Fortran.

## Functions
The following list describes the implemented functions, each validated to at least 8 digits of precision within the specified domains. Outside these domains, the same precision is not guaranteed.

* Exponential integral $\mathrm{Ei}(x)$
    * $\lbrace x \in \mathbb{R} \mid x \neq 0 \rbrace$
* Exponential integral $\mathrm{E}_1(x)$
    * $\lbrace x \in \mathbb{R} \mid x \neq 0 \rbrace$
* Exponential integral $\mathrm{E}_1(z)$
    * $z \in \mathbb{C} \setminus \left( \lbrace z \in \mathbb{C} \mid \Re(z) \lt 0,\thinspace 0 \lt |\Im(z)| \lt 10^{-6} \rbrace \cup \lbrace 0 \rbrace \right)$

The list above will be updated as new functions are added and tested. Next in line are:

* Bessel functions of the first kind $J_0(x)$ and $J_1(x)$
* Bessel functions of the second kind $Y_0(x)$ and $Y_1(x)$
* Struve functions $\mathbf{H}_0(x)$ and $\mathbf{H}_1(x)$
* Hypergeometric function ${}_2F_1(a, b; c; x)$

## Tests
Tests are conducted by comparing the ColSpecF results with those of [mpmath], an arbitrary-precision numerical library. These tests ensure at least 8 digits of precision within the specified domains.

Testing routines are built using [test-drive], a standard Fortran unit testing framework.

## Documentation
The [API documentation][Docs Website] for this library is generated using [FORD] and is deployed and hosted on [ReadTheDocs].

## References
Fortran code for evaluating special functions is sourced from the following websites:

* Association for Computing Machinery. 2012. [Collected Algorithms][calgo]
* Jason Blevins. 2004. [Alan Miller's Fortran Software][jblevins]
* Commonwealth Scientific and Industrial Research Organisation. 2004. [Software from Alan J. Miller][csiro]

ACM publications presenting the original algorithms are listed below:

1. Kathleen A. Paciorek. 1970. Algorithm 385: Exponential integral Ei(x). Commun. ACM 13, 7 (July 1970), 446–447. <https://doi.org/10.1145/362686.362696>
2. Donald E. Amos. 1990. Algorithms 683: a portable FORTRAN subroutine for exponential integrals of a complex argument. ACM Trans. Math. Softw. 16, 2 (June 1990), 178–182. <https://doi.org/10.1145/78928.78934>

## License
ColSpecF is distributed under two licenses based on code origin:

- Code comprising original contributions by Rodrigo Castro is licensed under the MIT License.
- Code adapted from Collected Algorithms (CALGO), published by the Association for Computing Machinery (ACM), is subject to the [ACM Software License Agreement][acmlic].

Users must comply with the applicable license for each portion of the code. See the [License][License File] for full details.

<!-- links -->
<!-- Badges -->
[Fortran Website]: https://fortran-lang.org/
[Fortran Badge]: https://img.shields.io/badge/Fortran-734f96?logo=fortran&style=flat
[Docs Website]: https://colspecf.readthedocs.io/
[Docs Badge]: https://img.shields.io/readthedocs/colspecf?color=blue
[Test Workflow]: https://github.com/rodpcastro/colspecf/actions/workflows/CI.yml
[Test Badge]: https://github.com/rodpcastro/colspecf/actions/workflows/CI.yml/badge.svg
[Coverage Website]: https://app.codecov.io/gh/rodpcastro/colspecf
[Coverage Badge]: https://codecov.io/github/rodpcastro/colspecf/badge.svg
[License File]: https://github.com/rodpcastro/colspecf/blob/main/LICENSE
[License Badge]: https://img.shields.io/badge/License-orange
<!-- Introduction -->
[Special Functions]: https://www.britannica.com/science/special-function
<!-- Tests -->
[mpmath]: https://mpmath.org/
[test-drive]: https://github.com/fortran-lang/test-drive
<!-- Documentation -->
[FORD]: https://forddocs.readthedocs.io/
[ReadTheDocs]: https://about.readthedocs.com/
<!-- References -->
[acm]: https://www.acm.org/
[calgo]: https://calgo.acm.org/
[jblevins]: https://jblevins.org/mirror/amiller/
[csiro]: https://wp.csiro.au/alanmiller/
<!-- License -->
[acmlic]: https://www.acm.org/publications/policies/software-copyright-notice
