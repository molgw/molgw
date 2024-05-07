# Instructions to compile HDF5

1. Download hdf5-1.12.2.tar.gz from the official website: https://www.hdfgroup.org/downloads/hdf5 (one has to create a free account)

2. Open the tarball:
```bash
tar -zxvf hdf5-1.12.2.tar.gz
cd hdf5-1.12.2
```

3. Configure with the `--enable-fortran` option:
```bash
./configure --enable-fortran --prefix=/opt/apps/hdf5-1.12.2
```

4. Compile:
```bash
make -j # possible warnings during the compilation
make test # all tests should be PASSED
make install
```

# Compile MOLGW with HDF5

In the my_machine.arch do the following:
1. Specify `HDF5_ROOT`. From the example above, it would be `/opt/apps/hdf5-1.12.2`
2. Specify the `-DHAVE_HDF5` flag
3. Possibly need to add `LDFLAGS += -Wl,-rpath,${HDF5_ROOT}/lib/` for dynamic linking
