# lagtraj Lagragian simulations trajectories

## Project components

Three steps:

1. Download a domain for a given date (small for Eulerian simulations, big for Lagrangian)

2. Produce trajectory

3. Extract forcing profiles along the trajectory

```bash
$> python -m lagtraj.domains.download [input_domain.yaml] [start date (yyyy-mm-dd)] [end date (yyyy-mm-dd)] [--file directories.yaml] [--overwrite]
$> python -m lagtraj.trajectory.create [input_trajectory.yaml]
$> python -m lagtraj.forcings.create [input_forcings.yaml]
```


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
