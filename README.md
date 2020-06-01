# lagtraj Lagragian simulations trajectories


## Producing a Lagrangian forcing

There are three steps to making forcing profiles with lagtraj:

1. Download a domain for a given date-range (small for Eulerian simulations,
big for Lagrangian)

2. Produce trajectory

3. Extract forcing profiles along the trajectory


## 0. Getting started

`lagtraj` stores everything (source data, definitions for how domains,
trajectories and forcings are set up) in a *data directory* (by default this
will be `data/` relative to where `lagtraj` is invoked). The directory
structure is as follows:

```bash
data
├── domains
│   ├── eurec4a_circle_eul.yaml
│   └── eurec4a_circle_eul_data
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

You can either make your own domain/forcing/trajectory definition (these are
stored in yaml-files) by creating a `meta.yaml` file in the relevant directory
or use the ones that `lagtraj` comes with by running the following command:

```bash
$> python -m lagraj.input_examples

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

    $> python -m lagtraj.domain.download lagtraj://eurec4a_20191209_12_eul 2020/01/01 2020/01/08
```

## 1. Making source data available

`lagtraj` is based around making all data required for interpolation,
integration and forcing calculation being available before trajectories are
integrated. This was done to reduce the number of data requests required to the
data storage backends (e.g. ECMWF), but does mean that *the expected extent
that a trajetory will reach must been known before performining a trajectory
integration*, otherwise `lagtraj` will issue a warning when the edge of the
available domain is reached.

Either create your own domain defition in `data/domains/<domain_name>/meta.yaml` and run

```bash
$> python -m lagtraj.domain.download <domain_name> [start date (yyyy-mm-dd)] [end date (yyyy-mm-dd)]
```

Or use one of the domain defitions included with `lagtraj` (e.g.
`eurec4a_north_atlantic`


```bash
$> python -m lagtraj.domain.download lagtraj://eurec4a_north_atlantic [start date (yyyy-mm-dd)] [end date (yyyy-mm-dd)]
```


## 2. Producing a trajectory

```bash
$> python -m lagtraj.trajectory.create <trajectory_name>
```

## 3. Producing forcing profiles

```bash
$> python -m lagtraj.forcing.create <forcing_name>
```

# Implementation details

Required utilies:

- download ECMWF data

- conversion of physical variables

- interpolation (do we need Steffen interpolation?)

- smoothing

- plotting (to some extent)

- generating HighTune formatted output netCDF files

- timing (can use tqdm for this)


## Algorithmic approach:

a) Download all data needed (across entire domain, single/model levels+).
   Currently this creates daily files, but we may aim for the minimum number of 
   requests from ECMWF (However, grib to netcdf conversion via the cds api
   can only deal with files of up to about 10GB).

b) Split these into daily files

Considerations:

- may need whole domain for plotting

- how many download requests can we make with ECMWF's data server? One big
  request seems better than many small ones
