# Holodeck2 C++

NOTE: build must be completed using correct version of python that will later use the module!
CMake build seems to be working correctly to generate shared library.
Using `setup.py` file (`python setup.py build_ext --inplace`) also seems to work.
NOTE: if the shared library is built against a different version of python, the module load will return an empty module (*but not raise any error!*).

## To-Do

### BUGS / URGENT

[ ] Implement scatter is MMBulge relation.  NOTE: the current hacky fix for avoiding mass-ratios above unity will need to be updated!
[ ] Implement non-GW hardening.

### General Improvements 

[ ] Add wrapper python-module that provides access to cpp module.  Use this to ensure module is loaded correctly (i.e. from same environment as library was built under).
[ ] Use `doctest` for C++ unit testing!
[ ] Change to better units (i.e. pc, Msol, Myr)
[ ] Handle low-mass MBHs going to un-physically high-densities.  Perhaps, implement occupancy fraction to zero the number-density of low-mass MBHs
[ ] Explore RNG Optimization.  The Poisson samples are extremely expensive (currently dominating runtime for GWB calculations; using boost::random::poisson_distribution).  It's possible that (at times?) a more custom approach could be better, for example, since each poisson_distribution is being samples `num_realizations` times, it might be faster to precompute all of the uniform random numbers, and then simultaneously invert the poisson CDF to find the corresponding samples for all of them together.  This would avoid searching the full CDF each time.  However, it looks like boost only uses inversion sampling for certain values of the mean, so this might only produce a speed up in those cases.  **First step**: do some simple tests about whether or not optimization is possible.

