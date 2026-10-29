#ifndef CURVATURE_H
#define CURVATURE_H

// Spin polarization at first order in the gradients of beta (and zeta) including the first
// non-vanishing correction in the derivatives of the normal vector (third order of the
// q-expansion, N = 3) of arXiv:2509.14301, trivial intersection ybar = x only.
// The geometry (derivatives of the fitted function F up to third order, and dSigma aligned
// with the analytic normal) is read from the binary file written by write_geometry.py.
// Formulas: curvature_correction.pdf; Python reference: spin_curvature.py.

#include <array>
#include <fstream>
#include <string>
#include <vector>

#include "particle_data_group.h"
#include "surface.h"

struct geo_cell {
    double Fa[4];             // d_a F             (Cartesian, lower indices)
    double Fab[4][4];         // d_a d_b F
    double Fabc[4][4][4];     // d_a d_b d_c F
};

// Reads the binary geometry file, checks that it belongs to `hypersurface` (same number of cells,
// same positions) and replaces every cell.dsigma with the dSigma aligned with the analytic normal.
void read_geometry(const std::string &filename, std::vector<element> &hypersurface,
                   std::vector<geo_cell> &geometry);

// Mean spin vector S^mu(p) for momentum (pT, phi, y_rap).  Cells with |p.n|/m < cut are left out of
// the N = 3 term (grazing world lines, where it is singular); eq. (46) and N_p always use all cells.
// Writes one line (28 columns):
//   pT phi y  N_p  S1v^0..3 S1s^0..3  S3v^0..3 S3s^0..3  S3litv^0..3 S3lits^0..3
//   N_p   : sum over cells of (p.dSigma) n_F, the same quantity as the "Denominator" column of the
//           improved_polarization_* outputs
//   S1v, S1s : eq. (46) (N = 1), thermal-vorticity and thermal-shear contributions to the mean spin
//           vector S^mu(p), ALREADY DIVIDED by N_p: S1v = P_vorticity/Denominator and
//           S1s = P_shear/Denominator of improved_polarization_* (dimensionless)
//   S3v, S3s : complete N = 3 correction (total derivatives through ybar(x), eq. (25)), vorticity and
//           shear contributions, divided by the same N_p
//   S3litv, S3lits : N = 3 term as written in eq. (47) (Delta(ybar) not differentiated), vorticity and
//           shear contributions, for comparison only.  Unlike S3v, S3litv does not vanish: the literal
//           form fails the Killing-field check (rigid rotation must give no correction)
// "Vorticity" and "shear" = contributions of the antisymmetric and symmetric parts of d_m beta_n;
// the formulas are linear in d beta, so the parts add up exactly.  Only the sum is free of
// gradients normal to the surface (see the discussion after eq. (46)).
// Mean spin vector including the correction: S1v + S1s + S3v + S3s  (polarization P = 2 S for
// spin 1/2).  Multiply by N_p to obtain numerators in the convention of the other primary files.
void curvature_polarization(double pT, double phi, double y_rap, pdg_particle particle,
                            std::vector<element> &hypersurface, std::vector<geo_cell> &geometry,
                            double cut, bool theta, std::ofstream &fileout);

#endif
