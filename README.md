# lagtraj Lagragian simulations trajectories

![lagtraj](https://github.com/EUREC4A-UK/lagtraj/workflows/lagtraj/badge.svg)

![trajectory example](docs/eurec4a_20191209_12_lag.png)


## Producing a Lagrangian forcing

There are four steps to making forcing profiles with lagtraj:

1. Download a domain for a given date-range (small for Eulerian simulations,
big for Lagrangian)

2. Produce trajectory

3. Extract forcing profiles along the trajectory

4. Convert forcing to desired output format

## 0. Getting started

### Installing lagtraj

`lagtraj` (and all its dependencies) can be installed with pip directly
from github:

```bash
$> python -m pip install git+https://github.com/EUREC4A-UK/lagtraj
```

`lagtraj` requires Python 3 and is tested with `python3.6` but later
versions should work too.

Once installed all `lagtraj`'s commands are available from any directory
and the follow the pattern

```bash
$> python -m lagtraj.<command>
```

### lagtraj input and output

`lagtraj` stores everything (source data, definitions for how domains,
trajectories and forcings are set up) in a *data directory* (by default this
will be `data/` relative to where `lagtraj` is invoked). The directory
structure is as follows:

```bash
data
├── domains
│   ├── eurec4a_circle.yaml
│   └── eurec4a_circle_data
│       ├── an_model_2020-01-01.nc
│       :
│       └── fc_single_2020-01-03.nc
├── forcings
│   ├── eure4a_20191209_12_eul.yaml
│   └── eure4a_20191209_12_eul.nc
└── trajectories
    ├── eure4a_20191209_12_eul.yaml
    └── eure4a_20191209_12_eul.nc
```

The `name` of each domain/trajectory/forcing inside `lagtraj` will be the
full filename without the `.yaml`-extension. E.g. the domain definition in
`domains/eurec4a_circle.yaml` will have the name `eurec4a_circle` inside
`lagtraj`.

You can either make your own domain/forcing/trajectory definition (these
are stored in yaml-files) by creating a yaml-file in the relevant
directory or use one that `lagtraj` comes with. You can list the
input-defintions bundled with your copy of `lagtraj` by running the
following command:

```bash
$> python -m lagraj.input_definitions.examples

The following domain/forcing/trajectory definitions are currently included
with lagtraj:

lagtraj://
 ├── domains
 │   ├── eurec4a_north_atlantic
 │   └── eurec4a_circle_eul
 ├── forcings
 │   └── eurec4a_20191209_12_eul
 └── trajectories
     ├── eurec4a_20191209_12_lin
     └── eurec4a_20191209_12_eul


To use for example the `eurec4a_north_atlantic` domain definition
for downloading domain data run lagtraj.domain.download as follows:

    $> python -m lagtraj.domain.download lagtraj://eurec4a_circle 2020/01/01 2020/01/08
```

## 1. Making source data available

`lagtraj` is based around making all data required for interpolation,
integration and forcing calculation available before trajectories and forcings
are calculated. This was done to reduce the number of data requests required
to the data storage backends (e.g. ECMWF), but does mean that *the expected
extent that a trajectory will reach must been known before performining
a trajectory integration*, otherwise `lagtraj` will issue a warning when the
edge of the available domain is reached.

Either create your own domain definition in `data/domains/<domain_name>.yaml` and run

```bash
$> python -m lagtraj.domain.download <domain_name> [start date (yyyy-mm-dd)] [end date (yyyy-mm-dd)]
```

Or use one of the domain definitions included with `lagtraj` (e.g.
`eurec4a_north_atlantic`


```bash
$> python -m lagtraj.domain.download lagtraj://eurec4a_circle [start date (yyyy-mm-dd)] [end date (yyyy-mm-dd)]
```

An optional flag `--retry-rate <num_minutes>` causes `lagtraj` to continue
retrying download of submitted data requests every `num_minutes` minutes until
all data has been downloaded. Every time this command is run it will attempt to
download only data not yet downloaded.


## 2. Producing a trajectory

Once you have downloaded the required domain data you can either create
a trajectory definition in `data/trajectories/<trajectory_name>.yaml` and run

```bash
$> python -m lagtraj.trajectory.create <trajectory_name>
```

Or use one of the domain definitions included with `lagtraj` (e.g.
`eurec4a_20191209_12_lag`


```bash
$> python -m lagtraj.trajectory.create lagtraj://eurec4a_20191209_12_lag
```

The created trajectory will be stored in `data/trajectories/<trajectory_name>.nc`.

## 3. Producing forcing profiles

```bash
$> python -m lagtraj.forcing.create <forcing_name>
```

## 4. Convert forcing to desired output format

```bash
$> python -m lagtraj.forcing.create <forcing_name> <conversion_specification>
```

e.g.

```bash
$>  python -m lagtraj.conversion.era5 eurec4a_20200202_12_lag lagtraj://eurec4a_dephy
$>  python -m lagtraj.conversion.era5 eurec4a_20200202_12_lag lagtraj://eurec4a_kpt
```

# Contributing and comments

Please feel free to [open an
issue](https://github.com/EUREC4A-UK/lagtraj/issues/new) if you have any
comments/questions/issues about `lagtraj`. Thank you!
