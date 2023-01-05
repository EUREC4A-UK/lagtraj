# lagtraj Lagragian simulations trajectories

![lagtraj](https://github.com/EUREC4A-UK/lagtraj/workflows/lagtraj/badge.svg)

![trajectory example](docs/eurec4a_20191209_12_lag.png)


## Producing a Lagrangian forcing

There are three steps to making forcing profiles with lagtraj:

1. [Download source data domain](#1-making-source-data-available) for a given
   date-range (small for Eulerian simulations, big for Lagrangian)

2. [Produce trajectory](#2-producing-a-trajectory)

3. [Extract forcing
   profiles](#3-producing-forcing-profiles)
   along the trajectory (with optional conversion to target a specific LES/GCM
   model)

The guide below first details how to install lagtraj and then guides your each
of these three steps.


## 0. Getting started

### Installing lagtraj

The most recent tagged version of `lagtraj` (and all its dependencies) can be
installed with pip directly from pipy:

```bash
$> python -m pip install lagtraj
```

*NOTE: if you are intending to modify `lagtraj` yourself you should check out
the [development notes](docs/developing.md).*

`lagtraj` requires Python 3 and is developed and tested with `python3.8` (in
that we aim to follow the recommendations of
[NEP29](https://numpy.org/neps/nep-0029-deprecation_policy.html)) but later
versions should work too (it may work with earlier versions too).

Once installed all `lagtraj`'s commands are available from any directory
and the follow the pattern

```bash
$> python -m lagtraj.<command>
```

### lagtraj input and output

`lagtraj` stores everything (both source data and *input definitions*
describing how domains, trajectories and forcings are set up) in a *data
directory* (by default this will be `data/` relative to where `lagtraj` is
invoked). The directory structure is as follows:

```bash
data
├── domains
│   ├── eurec4a_circle.yaml
│   └── eurec4a_circle_data
│       ├── an_model_2020-01-01.nc
│       :
│       └── fc_single_2020-01-03.nc
├── forcings
│   ├── eure4a_20200103_lag_testcase.yaml
│   └── eure4a_20200103_lag_testcase.nc
└── trajectories
    ├── eure4a_20200103_lag_testcase.yaml
    └── eure4a_20200103_lag_testcase.nc
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
$> python -m lagtraj.input_definitions.examples
```

Which will print

```bash
The following domain/forcing/trajectory definitions are currently included
with lagtraj:

lagtraj://
 ├── domains
 │   ├── drydown_cardington_local
 │   ├── eurec4a_circle
 │   └── eurec4a_north_atlantic
 ├── forcings
 │   ├── drydown_cardington_20200420_00_eul
 │   ├── eurec4a_campaign_eulerian
 │   ├── eurec4a_20200128_first
 │   ├── eurec4a_20200202_first_short_press
 │   ├── eurec4a_20200202_first_short
 │   ├── eurec4a_20200202_first
...
```

## 1. Making source data available

`lagtraj` is based around making all data required for interpolation,
integration and forcing calculation available before trajectories and forcings
are calculated. This was done to reduce the number of data requests required
to the data storage backends (e.g. ECMWF), but does mean that *the expected
spatial extent that a trajectory will reach must been known before performining
a trajectory integration*, otherwise `lagtraj` will issue a warning when the
edge of the available domain is reached.

In order to download the ERA5 input data for `lagtraj`, you need an account with
the Copernicus Data Store. You will also need to install the CDS api, see the [api-howto](https://cds.climate.copernicus.eu/api-how-to).

Either create your own domain definition in `data/domains/<domain_name>.yaml` and run

```bash
$> python -m lagtraj.domain.download <domain_name> <start_date> <end_date>
```

Or use one of the domain definitions included with `lagtraj` (e.g.
`eurec4a_circle`


```bash
$> python -m lagtraj.domain.download lagtraj://eurec4a_circle <start_date> <end_date>
```
the `<start_date>` and `<end_date>` should be formatted as `YYYY/MM/DD`, e.g. `2020/02/02` for the 2nd of February 2020.

An optional flag `--retry-rate <num_minutes>` causes `lagtraj` to continue
retrying download of submitted data requests every `num_minutes` minutes until
all data has been downloaded. Every time this command is run it will attempt to
download only data not yet downloaded.

You can monitor the status of your requests via the [CDS requests page](https://cds.climate.copernicus.eu).
Download times for model level data on the CDS can be somewhat variable.

## 2. Producing a trajectory

Once you have downloaded the required domain data you can either create
a trajectory definition in `data/trajectories/<trajectory_name>.yaml` and run

```bash
$> python -m lagtraj.trajectory.create <trajectory_name>
```

Or use one of the trajectory definitions included with `lagtraj` (e.g.
`eurec4a_20200202_first_short`


```bash
$> python -m lagtraj.trajectory.create lagtraj://eurec4a_20200202_first_short
```

The created trajectory will be stored in `data/trajectories/<trajectory_name>.nc`.

## 3. Producing forcing profiles

To produce forcings you need to create a forcing definition in
`data/forcings/<forcing_name>.yaml` and run

```bash
$> python -m lagtraj.forcings.create <forcing_name> [--conversion <conversion_name>]
```

Or use one of the forcing definitions included with `lagtraj` (e.g.
`eurec4a_20200202_first_short`)

```bash
$> python -m lagtraj.forcings.create lagtraj://eurec4a_20200202_first_short [--conversion <conversion_name>]
```

### Forcing profiles conversion (targeting a specific GCM/LES)

When creating forcings it might be desirable to target a specific LES
(Large-Eddy Simulation) model or GCM (Global Circulation Model) by
converting the forcings to a specific format and setting parameters
specific to the model being targeted. This can be achieved by using the
`--conversion` flag and providing a `conversion_name`. `lagtraj` currently
comes bundled with functionality to target the
[KPT](https://www.lmd.jussieu.fr/~mpllmd/dephy2_forcages_communs/KPT_documentation.pdf)
LES and
[dephy](https://docs.google.com/document/d/118xP04jB9HO7Y2LqWk3HZpZ9n3CFujgzimLI7Ug8vO4)
LES format.

Conversion parameters are defined in a yaml-files similarly to how domain,
trajectory and forcings definitions are stored, with one important difference:
conversion definition files are associated with a specific forcing definition
file (i.e. each forcing conversion definition points to only one specific
forcing definition). `lagtraj` comes bundled with definitions for how to do
forcing conversion with sensible defaults that you can modify for each forcing
you wish to create.

To set the parameters for a conversion identifed by the name `kpt` for
converting a forcing with name `forcing_name` you should create a file in
`data/forcings/<forcing_name>.<conversion_name>.yaml`. Running a conversion
will the convert `data/forcings/<forcing_name>.nc` and save to
`data/forcings/<forcing_name>.<conversion_name>.nc`.

```bash
$> python -m lagtraj.forcings.create <forcing_name> [--conversion <conversion_name>]
```

Instead of creating a conversion definition starting from an empty file you can
bootstrap the process by using the default parameters for a given target model
included with lagtraj. This is achieved by using the `lagtraj://`-prefix when
choosing the conversion name. E.g. to create the forcing named
`eurec4a_20200202_first_short` bundled with `lagtraj` and have it converted to
the `dephy` format with the default parameters you would run

```bash
$> python -m lagtraj.forcings.create lagtraj://eurec4a_20200202_first_short --conversion lagtraj://dephy
```

This will create the un-converted forcing in
`data/forcings/eurec4a_20200202_first_short.nc`, the converted one in
`data/forcings/eurec4a_20200202_first_short.dephy.nc` and the default conversion definition
for targeting `dephy` will be copied to
`data/forcings/eurec4a_20200202_first_short.dephy.yaml`. You can then modify
the forcing parameters (for example change the number of levels) by editing
`data/forcings/eurec4a_20200202_first_short.dephy.yaml` and rerunning the
forcing creation with your local copy of the conversion definition (note the
absence of the `lagtraj://` prefix):

```bash
$> python -m lagtraj.forcings.create lagtraj://eurec4a_20200202_first_short --conversion dephy
```

You are of course welcome to rename the conversion however you like if for
example you'd like to have multiple different converted version of the same
forcings file.

# Contributing and comments

If you have any comments/questions/issues about `lagtraj` please feel free to
[open an issue](https://github.com/EUREC4A-UK/lagtraj/issues/new) or have a
look in [docs/developing.md](docs/developing.md) where we have added some notes
on how to get started. All contributions are very welcome!
