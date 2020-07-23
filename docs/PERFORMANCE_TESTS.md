# Performance tests for the Lagtraj code

The test data required for a quick performance test (5 hour time integration) can de downloaded as follows:

```bash
python -m lagtraj.domain.download lagtraj://eurec4a_circle_eul 2020/02/02 2020/02/02
```

## Speed testing

This is done with pyinstrument:

[https://pypi.org/project/pyinstrument/0.10.1/](https://pypi.org/project/pyinstrument/0.10.1/)

Trajectory creation can be tested with the following command

```bash
pyinstrument -m lagtraj.trajectory.create lagtraj://eurec4a_20200203_12_lag
```

## Memory usage

This can be done with python memory profiler

[https://pypi.org/project/memory-profiler/](https://pypi.org/project/memory-profiler/)

```bash
mprof run python -m lagtraj.trajectory.create lagtraj://eurec4a_20200203_12_lag
mprof plot
```

Memory usage during trajectory creating for the test case has so far been about 175MB.
