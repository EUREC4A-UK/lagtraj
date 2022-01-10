# Changelog

## [v0.1.0](https://github.com/EUREC4A-UK/lagtraj/tree/v0.1.0)

[Full Changelog](https://github.com/EUREC4A-UK/lagtraj/compare/...v0.1.0)

First tagged version of lagtraj!

This version provides complete functionality to produce large-scale forcings for
running Lagrangian (flow-following) limited area simulations (e.g. Large-Eddy
Simulations or Single Column Models) to produce convective clouds based on
ECMWF ERA5 reanalysis data. Output forcings can currently target the models
supporting the
[KPT](https://www.lmd.jussieu.fr/~mpllmd/dephy2_forcages_communs/KPT_documentation.pdf)
(e.g. [DALES](https://github.com/dalesteam/dales)) and
[DEPHY](https://docs.google.com/document/d/118xP04jB9HO7Y2LqWk3HZpZ9n3CFujgzimLI7Ug8vO4/edit)
(e.g. [MONC](https://github.com/Leeds-MONC/monc)) input formats.

Features:

- calculation of large-scale forcings along a lat/lon trajectory and conversion
  to target specific LES models
- calculation of air-mass trajectories following either a) fixed-pressure or b)
  fixed-height
- automatic download of ECMWF ERA5 reanalysis data from CDSAPI (handling both
  data requests, download and consistency checks)
- rudimentary support for calculation of forcing profiles over land (by
  exclusion of fixed-height profile points which are below the surface
  elevation)
