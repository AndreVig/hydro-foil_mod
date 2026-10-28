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

`./foil beta.dat OUTDIR --curv beta.geo [--cut c] [-R] [--point pT phi y ...]`

Computes eq. (46) of arXiv:2509.14301 (trivial intersection ybar = x) **plus** the first non-vanishing correction in the
derivatives of the normal vector (third order of the q-expansion, N = 3), using the analytic description of the
hypersurface. `beta.geo` is written for this `beta.dat` by `write_geometry.py` (analytic-fit tools); it stores, per cell,
the derivatives of the fitted function F up to third order and dSigma aligned with the analytic normal.

Output `OUTDIR/primary_curvature_(mid)rapidity`, one line per momentum:
`pT phi y N_p S1^0..3 S3^0..3 S3lit^0..3`
- `S1`: eq. (46), identical to (P_vorticity + P_shear)/Denominator of `--impr-std` on the same dSigma;
- `S3`: complete N = 3 correction; the result with the correction is `S1 + S3`;
- `S3lit`: N = 3 term as written in eq. (47) (for comparison only).

`--cut c` drops cells with |p.n|/m < c from the N = 3 term, which is singular for world lines tangent to the surface;
quote results for several cuts and for two fits of the surface. Adding `--impr-theta` excludes cells with p.dSigma < 0
from everything, as in the `--impr-theta` mode (output file name gets `_theta`). Feed-down is not included in this mode.
Formulas: `curvature_correction.pdf`; Python reference: `spin_curvature.py`.
