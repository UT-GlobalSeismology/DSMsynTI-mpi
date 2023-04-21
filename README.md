**********  Welcome to the DSM world.  **********

This package is a bundle of 3 separate programs: tish, tipsv, and spcsac.
You must build each program seperately.
```
% cd tish-mpi
% make
% cd ../tipsv-mpi
% make
% cd ../spcsac
% make
```

# Important note
Some workers compute only the toroidal contribution to the transverse component, and the spheroidal contribution to the radial component. However, note that it is important not to omit the toroidal contribution to the radial displacement or the spheroidal contribution to the transverse component, or errors on the order of 1% can occur.

# Papers describing the methods used by this software
- Kawai, K., N. Takeuchi, and R.J. Geller, Complete synthetic seismograms up to 2 Hz for transversely isotropic spherically symmetric media, Geophys. J. Int., 164, 411-424, 2006.
- Takeuchi, N., R.J. Geller, and P.R. Cummins, Highly accurate P-SV complete synthetic seismograms using modified DSM operators, Geophys. Res. Lett., 23, 1175-1178, 1996.
- Geller, R.J., and N. Takeuchi, A new method for computing highly accurate DSM synthetic seismograms, Geophys. J. Int., 123, 449-470, 1995.
- Cummins, P.R., R.J. Geller, T. Hatori, and N. Takeuchi, DSM complete synthetic seismograms: SH, spherically symmetric, case, Geophys. Res. Lett., 21, 533-536, 1994.
- Cummins, P.R., R.J. Geller, and N. Takeuchi, DSM complete synthetic seismograms: P-SV, spherically symmetric, case, Geophys. Res. Lett., 21, 1663-1666, 1994.


