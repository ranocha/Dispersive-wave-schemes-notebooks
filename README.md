# A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3908803.svg)](https://doi.org/10.5281/zenodo.3908803)

This repository contains some code used in the article
```
@article{ranocha2021broad,
  title={A Broad Class of Conservative Numerical Methods for Dispersive
         Wave Equations},
  author={Ranocha, Hendrik and Mitsotakis, Dimitrios and Ketcheson, David I},
  journal={Communications in Computational Physics},
  year={2021},
  month={02},
  volume={29},
  number={4},
  pages={979--1029},
  publisher={Global Science Press},
  doi={10.4208/cicp.OA-2020-0119},
  eprint={2006.14802},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

> We develop general tools to construct fully-discrete, conservative numerical methods and apply them to several nonlinear dispersive wave equations: Benjamin-Bona-Mahony (BBM), Fornberg-Whitham, Camassa-Holm, Degasperis-Procesi, Holm-Hone, and the BBM-BBM system. These full discretizations conserve all linear invariants and one nonlinear invariant for each system. The spatial semidiscretizations are built using the unifying framework of summation by parts operators and include finite difference, spectral collocation, and both discontinuous and continuous finite element methods. Classical time integration schemes such as Runge-Kutta methods and linear multistep methods are applied and the recent relaxation technique is used to enforce temporal conservation.


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
