<!--
Goodness-of-fit test for autoregressive functional processes and specification test for OU process from a functional perspective. Software companion for "A goodness-of-fit test for functional time series: a specification test to SDE"
Authors: Alejandra López-Pérez and Javier Álvarez Liébana (@DadosDeLaplace)
-->

Goodness-of-fit test for autoregressive functional processes and specification test for OU process from a functional perspective
======

[![License:
GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<!-- <img src="" alt="goffda  hexlogo" align="right" width="200" style="padding: 0 15px; float: right;"/> -->

Overview
-------

Software companion for the article [A goodness-of-fit test for functional time series: a specification test to SDE (López-Perez, Álvarez-Liébana, González-Manteiga and Febrero Bande, 2021)](https://arxiv.org/), currently submitted in JOURNAL. It implements the proposed estimators and goodness-of-fit tests for 
autoregressive functional processes of order one in Hilbert spaces (ARH(1) processes), as well as the specification test for diffusion processes, focused on testing [Ornstein-Uhlenbeck processes](https://www.sciencedirect.com/science/article/abs/pii/S016771521630044X).

It also allows to replicate simulations and the data application presented.

Installation
------------

Get the released version from GitHub:

``` r
# Install the package (devtools package is required to be installed)
devtools::install_github()

# Load package
library(goffda)
```


**Feel free to use** any of the contents but don't forget to cite them.

References
----------

García-Portugués, E., Álvarez-Liébana, J., Álvarez-Pérez, G. and
González-Manteiga, W. (2021). A goodness-of-fit test for the functional
linear model with functional response. *Scand J Statist.* 2021; 48: 502– 528.
<a href="https://doi.org/10.1111/sjos.12486" class="uri">https://doi.org/10.1111/sjos.12486</a>

López-Perez, A., Álvarez-Liébana, J., González-Manteiga, W. and Febrero-Bande, M. (2021). 
A goodness-of-fit test for functional time series: a specification test to SDE. *arXiv: *
<a href="https://arxiv.org" class="uri">https://arxiv.org</a>
