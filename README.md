**_Welcome to the DSM world._**

This is a software for computing synthetic seismograms in a spherically symmetric, transversely isotropic (TI) model using the Direct Solution Method (DSM).

This package is a bundle of 3 separate programs: tish, tipsv, and spcsac.
You must build each program separately.

```
% cd tish-mpi
% make
% cd ../tipsv-mpi
% make
% cd ../spcsac
% make
```

# Important note

If you compute only the toroidal contribution to the transverse component and the spheroidal contribution to the radial component, this will lead to errors on the order of 1% or more. Therefore you **must** include the toroidal contribution to the radial displacement and the spheroidal contribution to the transverse component.

# Papers describing the methods and theory used by this software

-   Kawai, K., N. Takeuchi, and R.J. Geller, Complete synthetic seismograms up to 2 Hz for transversely isotropic spherically symmetric media, Geophys. J. Int., 164, 411-424, 2006.
-   Takeuchi, N., R.J. Geller, and P.R. Cummins, Highly accurate P-SV complete synthetic seismograms using modified DSM operators, Geophys. Res. Lett., 23, 1175-1178, 1996.
-   Geller, R.J., and N. Takeuchi, A new method for computing highly accurate DSM synthetic seismograms, Geophys. J. Int., 123, 449-470, 1995.
-   Cummins, P.R., R.J. Geller, T. Hatori, and N. Takeuchi, DSM complete synthetic seismograms: SH, spherically symmetric, case, Geophys. Res. Lett., 21, 533-536, 1994.
-   Cummins, P.R., R.J. Geller, and N. Takeuchi, DSM complete synthetic seismograms: P-SV, spherically symmetric, case, Geophys. Res. Lett., 21, 1663-1666, 1994.
-   Geller, R.J., and T. Ohminato, Computation of synthetic seismograms and their partial derivatives for heterogeneous media with arbitrary natural boundary conditions using the Direct Solution Method, Geophys. J. Int., 116, 421-446, 1994.

# Authorship and copyright of software

This software was written and improved, and is copyrighted &copy;, by the members of the Global Seismology Group of the University of Tokyo from 1994 to present.

# License

This software is made available under the GNU Public License v3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
