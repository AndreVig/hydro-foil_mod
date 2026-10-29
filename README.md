# HYDROdynamic Freeze-Out IntegraL - hydro-foil
Compute integrals on the freeze-out hypersurface. By default, the integration happens parallelly on NTHREADS = number_of_threads_available - 2. This number can be modified prior to compiling in "integrals.cpp". The use of OpenMP can be completely switched off by commenting the "OPENMP" flag in the makefile

 The typical execution from command line is:
./foil PATH_TO/beta.dat PATH_TO/outputfile

The main file should be modified changing the function to use depending on the calculation. The functions available at the moment are:

* **void polarization_midrapidity(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);** Computes the mean spin vector for particles emitted at midrapidity for the particle "particle" using the freezeout file "freezeout_sup". For the thermal vorticity, eq. 64 of 2304.02276 is used. Refer also to 2103.10917 for details about the formulae.
The output is written in "fileout" as a table:

| pt | phi | denominator | numerator_varpi | numerator_xi |

* **void polarization_midrapidity_linear(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);** Same as the previous function but uses the linear approximation for the vorticity-induced polarization. This function is faster.



* **void polarization_projected(double pT, double phi, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);** Similar to "polarization_midrapidity_linear" but computes polarization projecting the gradients tangent to the hypersurface at point x before integrating.  

* **void polarization_rapidity(double pT, double phi, double y_rap, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);** Same as "polarization_midrapidity", but the table now includes the rapidity "y":

| pt | phi | y | denominator | numerator_varpi | numerator_xi |


* **void spectrum_rapidity(double pT, double phi, double y_rap, pdg_particle particle, vector<element> &freeze_out_sup, ofstream &fileout);** Integrates the thermal spectrum of "particle" using the Cooper-Freye prescription.



## Curvature correction (--curv)

`./foil beta.dat OUTDIR --curv beta.geo [--cut c] [-R] [--point pT phi y ...] [--pt-grid n min max] [--phi-grid n] [--y-grid n min max]`

Computes eq. (46) of arXiv:2509.14301 (trivial intersection ybar = x) **plus** the first non-vanishing correction in the
derivatives of the normal vector (third order of the q-expansion, N = 3), using the analytic description of the
hypersurface. `beta.geo` is written for this `beta.dat` by `write_geometry.py` (analytic-fit tools); it stores, per cell,
the derivatives of the fitted function F up to third order and dSigma aligned with the analytic normal.

Output `OUTDIR/primary_curvature_(mid)rapidity`, one line per momentum (28 columns):
`pT phi y N_p S1v^0..3 S1s^0..3 S3v^0..3 S3s^0..3 S3litv^0..3 S3lits^0..3`
- `N_p`: the same quantity as the Denominator column of the other modes;
- `S1v`, `S1s`: eq. (46), thermal-vorticity and thermal-shear contributions to the mean spin vector, already divided by
  N_p, i.e. identical to P_vorticity/Denominator and P_shear/Denominator of `--impr-std` on the same dSigma;
- `S3v`, `S3s`: the N = 3 correction, split in the same way and divided by the same N_p;
- `S3litv`, `S3lits`: N = 3 term as written in eq. (47), split in the same way (for comparison only).
  `S3v` vanishes identically (a Killing field, i.e. rigid rotation, gives no correction); `S3litv` does not, which shows
  that the literal eq. (47) misses terms.

The mean spin vector with the correction is `S1v + S1s + S3v + S3s` (multiply by N_p for numerators). The split is by the
antisymmetric and symmetric parts of d_m beta_n; only the sum is free of gradients normal to the surface.

`--cut c` drops cells with |p.n|/m < c from the N = 3 term, which is singular for world lines tangent to the surface;
quote results for several cuts and for two fits of the surface. Adding `--impr-theta` excludes cells with p.dSigma < 0
from everything, as in the `--impr-theta` mode (output file name gets `_theta`). Feed-down is not included in this mode.
Momentum grid (default as in the other modes: 30 pT in [0, 6] GeV x 30 phi in [0, 2pi] x 20 y in [-1, 1] with `-R`):
`--pt-grid n min max` (n values, ends included), `--phi-grid n` (n uniform values 0 ... 2pi(n-1)/n, no duplicate of 0),
`--y-grid n min max` (with `-R`). The N = 3 term is not smooth in phi (cells crossing the cut), so do not go below
about 24 phi values when S3 matters.
Formulas: `curvature_correction.pdf`; Python reference: `spin_curvature.py`.
