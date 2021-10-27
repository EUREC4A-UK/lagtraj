# lagtraj Lagragian simulations trajectories

![lagtraj](https://github.com/EUREC4A-UK/lagtraj/workflows/lagtraj/badge.svg)

![trajectory example](docs/eurec4a_20191209_12_lag.png)


## Producing a Lagrangian forcing

There are four steps to making forcing profiles with lagtraj:

1. Download a domain for a given date-range (small for Eulerian simulations,
big for Lagrangian)

2. Produce trajectory

3. Extract forcing profiles along the trajectory (with optional conversion to
   target a specific LES/GCM model)

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

    $> python -m lagtraj.domain.download lagtraj://eurec4a_20191209_12_eul 2020/01/01 2020/01/08
```

## 1. Making source data available

`lagtraj` is based around making all data required for interpolation,
integration and forcing calculation being available before trajectories are
integrated. This was done to reduce the number of data requests required to the
data storage backends (e.g. ECMWF), but does mean that *the expected extent
that a trajectory will reach must been known before performining a trajectory
integration*, otherwise `lagtraj` will issue a warning when the edge of the
available domain is reached.

Either create your own domain definition in `data/domains/<domain_name>.yaml` and run

```bash
$> python -m lagtraj.domain.download <domain_name> [start date (yyyy-mm-dd)] [end date (yyyy-mm-dd)]
```

Or use one of the domain definitions included with `lagtraj` (e.g.
`eurec4a_north_atlantic`


```bash
$> python -m lagtraj.domain.download lagtraj://eurec4a_north_atlantic [start date (yyyy-mm-dd)] [end date (yyyy-mm-dd)]
```


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

## 3. Producing forcing profiles

To produce forcings you need to create a forcing definition in
`data/forcings/<forcing_name>.yaml` and run

```bash
$> python -m lagtraj.forcing.create <forcing_name> [--conversion <conversion_name>]
```

Or use one of the forcing definitions included with `lagtraj` (e.g.
`eurec4a_20200202_12_lag`)

```bash
$> python -m lagtraj.forcing.create lagtraj://eurec4a_20200202_12_lag [--conversion <conversion_name>]
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
file. To set the parameters for a conversion identifed by the name `kpt` for
converting a forcing with name `forcing_name` you should create a file in
`data/forcings/<forcing_name>.<conversion_name>.yaml`. Running a conversion
will the convert `data/forcings/<forcing_name>,nc` and save to
`data/forcings/<forcing_name>.<conversion_name>.nc`.

```bash
$> python -m lagtraj.forcing.create <forcing_name> [--conversion <conversion_name>]
```

Instead of creating a conversion definition starting from an empty file you can
bootstrap the process by using the default parameters for a given target model
included with lagtraj. This is achieved by using the `lagtraj://`-prefix when
choosing the conversion name. E.g. to create the `eurec4a_20200202_12_lag`
bundled with `lagtraj` and have it converted to the `dephy` format with the
default parameters you would run

```bash
$> python -m lagtraj.forcing.create lagtraj://eurec4a_20200202_12_lag --conversion lagtraj://dephy
```

This will create the un-converted forcing in
`data/forcings/eurec4a_20200202_12_lag.nc`, the converted on in
`data/forcings/eurec4a_20200202_12_lag.dephy.nc` and the conversion definition
will be copied to `data/forcings/eurec4a_20200202_12_lag.dephy.yaml`. You can
then modify the forcing parameters (for example change the number of levels) by
editing `data/forcings/eurec4a_20200202_12_lag.dephy.yaml` and rerunning the
forcing creation with your local copy of the conversion definition (note the
absence of the `lagtraj://` prefix):

```bash
$> python -m lagtraj.forcing.create lagtraj://eurec4a_20200202_12_lag --conversion dephy
```

You are of course welcome to rename the conversion however you like if for
example you'd like to have multiple different converted version of the same
forcings file.


## Modifying lagtraj

If you spot a bug in lagtraj or would like to add features yourself please have
a look in [docs/developing.md](docs/developing.md) where we have added some
notes on how to get started. All contributions are very welcome!
