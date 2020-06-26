# A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3908803.svg)](https://doi.org/10.5281/zenodo.3908803)

This repository contains some code used in the article
```
@online{ranocha2020broad,
  title={A Broad Class of Conservative Numerical Methods for Dispersive
         Wave Equations},
  author={Ranocha, Hendrik and Mitsotakis, Dimitrios and Ketcheson, David I},
  year={2020},
  month={06},
  eprint={2006.XXXXX},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

> We further develop general tools to construct discretely conservative numerical methods and apply them to several dispersive wave equations: Benjamin-Bona-Mahony (BBM), Fornberg-Whitham, Camassa-Holm, Degasperis-Procesi, Holm-Hone, and the BBM-BBM system. These full discretizations conserve all linear invariants and one nonlinear invariant for each system. The spatial semidiscretizations are built using the unifying framework of summation by parts operators and include finite difference, spectral collocation, and both discontinuous and continuous finite element methods. Classical time integration schemes such as Runge-Kutta methods and linear multistep methods are applied and the recent relaxation technique is used to enforce temporal conservation.


If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please cite this repository as
```
@misc{ranocha2020broadRepro,
  title={{Dispersive-wave-schemes-notebooks}.
         {A} Broad Class of Conservative Numerical Methods for Dispersive
         Wave Equations},
  author={Ranocha, Hendrik and Mitsotakis, Dimitrios and Ketcheson, David I},
  year={2020},
  month={06},
  howpublished={\url{https://github.com/ranocha/Dispersive-wave-schemes-notebooks}},
  doi={10.5281/zenodo.3908803}
}
```

The numerical experiments were carried out on a workstation running Kubuntu 20.04. The MATLAB scripts were executed using MATLAB R2018a. The Jupyter notebooks running Julia kernels used Julia v1.4.2.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
