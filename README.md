<p align="center">
    <img src="https://raw.githubusercontent.com/rodpcastro/colspecf/refs/heads/main/media/colspecf_logo.svg", alt="ColSpecF Logo">
</p>

[![Fortran][Fortran Badge]][Fortran Website]
[![Docs][Docs Badge]][Docs Website]
[![Test][Test Badge]][Test Workflow]
[![Coverage][Coverage Badge]][Coverage Website]

# ColSpecF
ColSpecF (Collected Special Functions) is a [Fortran][Fortran Website] library for evaluating mathematical [Special Functions], built around adaptations of algorithms collected from [several sources](#references).

## Functions
The following list describes the implemented functions and their domains.

* Exponential integrals $\mathrm{Ei}(x)$ and $\mathrm{E}_1(x)$
    * $\lbrace x \in \mathbb{R} \mid x \neq 0 \rbrace$
* Exponential integral $\mathrm{E}_1(z)$
    * $\lbrace z \in \mathbb{C} \mid x \neq 0 \rbrace$
* Bessel functions of the first kind $J_0(x)$ and $J_1(x)$
    * $x \in \mathbb{R}$
* Bessel functions of the second kind $Y_0(x)$ and $Y_1(x)$
    * $\lbrace x \in \mathbb{R} \mid x \gt 0 \rbrace$
* Gauss hypergeometric function ${}_2F_1(a, b; c; z)$
    * $a,\thinspace b,\thinspace c,\thinspace z \in \mathbb{C}$

The list above will be updated as new functions are added and tested. Next in line are:

* Struve functions $\mathbf{H}_0(x)$ and $\mathbf{H}_1(x)$

## Tests
Tests are conducted by comparing the ColSpecF results with those of [mpmath], an arbitrary-precision numerical library.

Testing routines are built using [test-drive], a standard Fortran unit testing framework.

Test results can be found [here][test_results].

## Documentation
The [API documentation][Docs Website] for this library is generated using [FORD] and is deployed and hosted on [ReadTheDocs].

## References
Fortran code for evaluating special functions is sourced from the following websites:

* Association for Computing Machinery. 2012. [Collected Algorithms][calgo]
* Jason Blevins. 2004. [Alan Miller's Fortran Software][jblevins]
* Commonwealth Scientific and Industrial Research Organisation. 2004. [Software from Alan J. Miller][csiro]
* Elsevier. 2025. [Elsevier Data Repository][elsvdata]
* John Burkardt. 2025. [Fortran77 Source Codes][jbf77]
* John Burkardt. 2025. [Fortran90 Codes][jbf90]

Publications presenting the original algorithms are listed below:

1. Kathleen A. Paciorek. 1970. Algorithm 385: Exponential integral Ei(x). Commun. ACM 13, 7 (July 1970), 446–447. <https://doi.org/10.1145/362686.362696>
2. Donald E. Amos. 1990. Algorithms 683: a portable FORTRAN subroutine for exponential integrals of a complex argument. ACM Trans. Math. Softw. 16, 2 (June 1990), 178–182. <https://doi.org/10.1145/78928.78934>
3. W. J. Cody. 1993. Algorithm 715: SPECFUN–a portable FORTRAN package of special function routines and test drivers. ACM Trans. Math. Softw. 19, 1 (March 1993), 22–30. <https://doi.org/10.1145/151271.151273>
4. N. Michel and M. V. Stoitsov. 2008. Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Pöschl-Teller-Ginocchio potential wave functions. Computer Physics Communications 178, 7 (April 2008), 535–551. <https://doi.org/10.1016/J.CPC.2007.11.007>

## License
ColSpecF is a Fortran library distributed under multiple licenses or permissions based on code origin. Users must comply with the applicable license or permission for each portion of the code. See the [License][License File] for full details.

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
<!-- Introduction -->
[Special Functions]: https://www.britannica.com/science/special-function
<!-- Tests -->
[mpmath]: https://mpmath.org/
[test-drive]: https://github.com/fortran-lang/test-drive
[test_results]: https://github.com/rodpcastro/colspecf/blob/gauss/test/test_results.md
<!-- [test_results]: https://github.com/rodpcastro/colspecf/blob/main/test/test_results.md -->
<!-- Documentation -->
[FORD]: https://forddocs.readthedocs.io/
[ReadTheDocs]: https://about.readthedocs.com/
<!-- References -->
[calgo]: https://calgo.acm.org/
[jblevins]: https://jblevins.org/mirror/amiller/
[csiro]: https://wp.csiro.au/alanmiller/
[elsvdata]: https://elsevier.digitalcommonsdata.com/
[jbf77]: https://people.sc.fsu.edu/~jburkardt/f77_src/f77_src.html
[jbf90]: https://people.sc.fsu.edu/~jburkardt/f_src/f_src.html
