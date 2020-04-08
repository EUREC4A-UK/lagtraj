# lagtraj Lagragian simulations trajectories

## Project components

Two steps:

1. Produce trajectory (could then be specified by other means)

2. Extract forcing profiles along trajectory

```bash
$> python -m lagtraj.produce.forcing_profiles [trajectory.nc]
$> python -m lagtraj.produce.lagrangian_trajectory -lat0 -lon0 -t_start
```


Required utilies:

- download ECMWF data

- conversion of physical variables

- interpolation (do we need Steffen interpolation?)

- smoothing

- plotting?

- generating HighTune formatted output netCDF files

- timing (can use tqdm for this)


## Algorithmic approach:

a) Download all data needed (across entire domain) at once in a single request
   to ECMWF

b) Or download interatively by each integration step, only downloading data
   required to do the next integration

Considerations:

- may need whole domain for plotting

- how many download requests can we make with ECMWF's data server? One big
  request many small ones
