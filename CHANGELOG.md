# Changelog

## [unreleased](https://github.com/EUREC4A-UK/lagtraj/tree/HEAD)

[Full Changelog](https://github.com/EUREC4A-UK/lagtraj/compare/v0.1.2...HEAD)

*improvements*

- Provide information on CDS and its API in the README file
  [\#194](https://github.com/EUREC4A-UK/lagtraj/pull/194) @sjboeing

*maintenance*

- add Zenodo integration and add citation file
  [\#195](https://github.com/EUREC4A-UK/lagtraj/pull/195) @leifdenby

- update versions of flake8, black and isort in pre-commit config and resolve
  linting issues arising from these updates
  [\#206](https://github.com/EUREC4A-UK/lagtraj/pull/206), @leifdenby

*features*
- add support for converting forcing files to target the SAM LES model
  (using the [SCAM IOP](https://www.cesm.ucar.edu/models/simple/scam) format)
  [\#160](https://github.com/EUREC4A-UK/lagtraj/pull/160) @xychen-ocn @pblossey

## [v0.1.2](https://github.com/EUREC4A-UK/lagtraj/tree/v0.1.2)

[Full Changelog](https://github.com/EUREC4A-UK/lagtraj/compare/v0.1.1...v0.1.2)

*bugfixes*

- fix text shown when using CLI command to list lagtraj bundled examples
  [\#187](https://github.com/EUREC4A-UK/lagtraj/pull/187) @leifdenby

- fix for forcing conversion input definition filepath (to remove
  `lagtraj://`-prefix) [\#185](https://github.com/EUREC4A-UK/lagtraj/pull/185)
  @leifdenby

- remove inaccurate description of how skin temperature is calculated in KPT
  formatted output files
  [\#182](https://github.com/EUREC4A-UK/lagtraj/pull/182/) @sjboeing

*maintenance*

- handle time-coordinate out-of-order domain data being returned from CDS for
  ERA5 domain data [\#189](https://github.com/EUREC4A-UK/lagtraj/pull/189)
  @leifdenby

- switch to using python v3.8 for continuous integration
  [\#186](https://github.com/EUREC4A-UK/lagtraj/pull/186) @leifdenby


## [v0.1.1](https://github.com/EUREC4A-UK/lagtraj/tree/v0.1.1)

[Full Changelog](https://github.com/EUREC4A-UK/lagtraj/compare/v0.1.0...v0.1.1)

*bugfixes*

- ensure package data (example input definitions and ERA5 level definitions)
  are installed when installing `lagtraj` from pypi
  [\#178](https://github.com/EUREC4A-UK/lagtraj/pull/178) @leifdenby

- fix bug preventing creation of KPT files related to change in ERA5 units for
  `sdor` variable
  [\#174](https://github.com/EUREC4A-UK/lagtraj/pull/174) @sjboeing & @leifdenby


## [v0.1.0](https://github.com/EUREC4A-UK/lagtraj/tree/v0.1.0)

[Full Changelog](https://github.com/EUREC4A-UK/lagtraj/compare/...v0.1.0)

First tagged version of lagtraj!

This version provides complete functionality to produce large-scale forcings for
running Lagrangian (flow-following) limited area simulations (e.g. Large-Eddy
Simulations or Single Column Models) to simulate convective clouds based on
ECMWF ERA5 reanalysis data. Output forcings can currently target the models
supporting the
[KPT](https://www.lmd.jussieu.fr/~mpllmd/dephy2_forcages_communs/KPT_documentation.pdf)
(e.g. [DALES](https://github.com/dalesteam/dales)) and
[DEPHY](https://docs.google.com/document/d/118xP04jB9HO7Y2LqWk3HZpZ9n3CFujgzimLI7Ug8vO4/edit)
(e.g. [MONC](https://github.com/Leeds-MONC/monc), DEPHY functionality yet on
MOSRS trunk) input formats.

Features:

- calculation of large-scale forcings along a lat/lon trajectory and conversion
  to target specific LES models
- calculation of air-mass trajectories following either a) fixed-pressure, b)
  fixed-height or fixed velocity from a starting point
- automatic download of ECMWF ERA5 reanalysis data from CDSAPI (handling both
  data requests, download and consistency checks)
- rudimentary support for calculation of forcing profiles over land (by
  exclusion of fixed-height profile points which are below the surface
  elevation)
