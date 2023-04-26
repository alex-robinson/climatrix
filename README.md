# climatrix
Library to provide transient climatic forcing for an ice-sheet model given several snapshots (interpolate over a matrix of geometry and anomalies)

## Quick start

The test program `tests/test_climatrix.f90` demonstrates the basic functionality of the climatrix library.
Hopefully it is self-explanatory. To get it running, perform the following steps.

### Configuration

Define a compiler configuration file in the directory `config/`. Essentially, provide the compiler name and the location of the installed NetCDF libraries for your system. Use `config/pagos_gfortran` as a template.Then, configure a Makefile for your system:

```python
python config.py config/pagos_gfortran
```

where `pagos_gfortran` is the name of your config file.

### Compilation

Compile the test program:

```
make clean
make test
```

This will compile the test program `tests/test_climatrix.f90` and store the executable as `libclimatrix/bin/test_climatrix.x`. 

### Run the program

Make sure you have the necessary input data stored in the `input/` directory. Then run the program:

```
./libclimatrix/bin/test_climatrix.x
```

That's it! 
