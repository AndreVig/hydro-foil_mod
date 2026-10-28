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
// Writes one line:  pT phi y  N_p  S1^0..3  S3^0..3  S3lit^0..3
//   S1    : eq. (46) (N = 1), same normalisation as (P_vorticity + P_shear)/Denominator of
//           improved_polarization_* (dimensionless)
//   S3    : complete N = 3 correction (total derivatives through ybar(x), eq. (25))
//   S3lit : N = 3 term as written in eq. (47) (Delta(ybar) not differentiated), for comparison
// The polarization including the correction is S1 + S3.
// theta = true: cells with p.dSigma < 0 are excluded from everything (as --impr-theta).
void curvature_polarization(double pT, double phi, double y_rap, pdg_particle particle,
                            std::vector<element> &hypersurface, std::vector<geo_cell> &geometry,
                            double cut, bool theta, std::ofstream &fileout);

#endif
