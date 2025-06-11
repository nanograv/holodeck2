# Holodeck2 C++

NOTE: build must be completed using correct version of python that will later use the module!
CMake build seems to be working correctly to generate shared library.
Using `setup.py` file (`python setup.py build_ext --inplace`) also seems to work.
NOTE: if the shared library is built against a different version of python, the module load will return an empty module (*but not raise any error!*).


## To-Do
[ ] Add wrapper python-module that provides access to cpp module.  Use this to ensure module is loaded correctly (i.e. from same environment as library was built under).
[ ] Use `doctest` for C++ unit testing!
[ ] Change to better units (i.e. pc, Msol, Myr)
[ ] Implement occupancy fraction to zero the number-density of low-mass MBHs



