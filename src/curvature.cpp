// Spin polarization with the first non-vanishing correction in the derivatives of the normal
// vector (N = 3 term), arXiv:2509.14301, ybar = x.  See curvature.h and curvature_correction.pdf.
//
// Notation (all Cartesian, metric diag(+,-,-,-), eps^{0123} = +1):
//   N_a = d_a F (any normalisation), Gamma_nu = N_nu/|p.N|, Delta^{ae} = g^{ae} - N^a p^e/(p.N)
//   Y(x) = x + t(x) p: crossing of the world line with the surface;  t_ab, t_abc its derivatives
//   A^mu_N = (-i)^(N+1)/N! eps^{mu rho lam tau} p_lam T^(N)_{a..tau} V_rho^nu_sigma E^{a..sigma}_nu
// The mass-shell integral of the delta' term is done exactly: p^0 is a dual number.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>

#ifdef OPEN_MP
#include <omp.h>
#endif

#include "curvature.h"
#include "utils.h"

using namespace std;

// ------------------------------------------------------------------ dual numbers (value, d/dp0)
namespace {

struct Dual {
    double v, d;
    Dual(double v_ = 0.0, double d_ = 0.0) : v(v_), d(d_) {}
};
inline Dual operator+(Dual a, Dual b) { return {a.v + b.v, a.d + b.d}; }
inline Dual operator-(Dual a, Dual b) { return {a.v - b.v, a.d - b.d}; }
inline Dual operator-(Dual a) { return {-a.v, -a.d}; }
inline Dual operator*(Dual a, Dual b) { return {a.v * b.v, a.v * b.d + a.d * b.v}; }
inline Dual operator*(double s, Dual a) { return {s * a.v, s * a.d}; }
inline Dual operator*(Dual a, double s) { return {s * a.v, s * a.d}; }
inline Dual operator/(Dual a, Dual b) { return {a.v / b.v, (a.d * b.v - a.v * b.d) / (b.v * b.v)}; }
inline Dual operator/(Dual a, double s) { return {a.v / s, a.d / s}; }
inline Dual &operator+=(Dual &a, Dual b) { a.v += b.v; a.d += b.d; return a; }

// Fermi-Dirac 1/(exp(E)+1), numerically stable, with its derivative
inline Dual fermi(Dual E) {
    double f = (E.v > 0) ? exp(-E.v) / (1.0 + exp(-E.v)) : 1.0 / (exp(E.v) + 1.0);
    return {f, -f * (1.0 - f) * E.d};
}

const double Gm[4] = {1.0, -1.0, -1.0, -1.0};
const int PAIRS[10][2] = {{0,0},{0,1},{0,2},{0,3},{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}};

struct CellResult {
    Dual pdS, nF;
    Dual A1[4], A3d[4], A3dp[4], A3d_lit[4], A3dp_lit[4];
    double grazing;   // |p.N| / (m sqrt|N.N|)
};

// eps^{mu rho lam tau} p_lam X_{rho tau}
inline void eps_contract(const Dual pl[4], const Dual X[4][4], Dual out[4]) {
    for (int m = 0; m < 4; m++) out[m] = Dual();
    for (auto e : non_zero_levi()) {
        out[e[0]] += double(e[4]) * pl[e[2]] * X[e[1]][e[3]];
    }
}

// all per-cell quantities for one 4-momentum p (p^0 dual)
void cell_terms(const Dual p[4], const element &c, const geo_cell &g, double mass, double mutot,
                CellResult &R) {
    Dual pl[4];
    for (int a = 0; a < 4; a++) pl[a] = Gm[a] * p[a];
    // measure, distribution
    R.pdS = Dual();
    for (int a = 0; a < 4; a++) R.pdS += c.dsigma[a] * p[a];
    double beta[4];
    for (int a = 0; a < 4; a++) beta[a] = Gm[a] * c.u[a] / c.T;          // beta_a (lower)
    Dual Earg = Dual(-mutot / c.T);
    for (int a = 0; a < 4; a++) Earg += beta[a] * p[a];
    R.nF = fermi(Earg);
    Dual one(1.0);
    Dual a0 = -(R.nF * (one - R.nF));
    Dual a2 = -(R.nF * (one - R.nF) * (one - 6.0 * R.nF + 6.0 * R.nF * R.nF)) / 24.0;

    const double *n = g.Fa;                  // N_a
    Dual pn;
    for (int a = 0; a < 4; a++) pn += n[a] * p[a];
    double s = (pn.v > 0) - (pn.v < 0);
    double nn = 0;
    for (int a = 0; a < 4; a++) nn += Gm[a] * n[a] * n[a];
    R.grazing = fabs(pn.v) / (mass * sqrt(fabs(nn)));
    Dual ipn = one / pn;

    // D[a][b] = d_a Y^b = delta_a^b - N_a p^b/(p.N);  Dup[a][e] = Delta^{ae}
    Dual D[4][4], Dup[4][4];
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++) {
            D[a][b] = Dual(a == b ? 1.0 : 0.0) - n[a] * p[b] * ipn;
            Dup[a][b] = Dual(a == b ? Gm[a] : 0.0) - (Gm[a] * n[a]) * p[b] * ipn;
        }
    // t_ab
    Dual t2[4][4];
    {
        Dual tmp[4][4];  // tmp[a][c] = sum_b D[a][b] dn[c][b]
        for (int a = 0; a < 4; a++)
            for (int cc = 0; cc < 4; cc++) {
                Dual acc;
                for (int b = 0; b < 4; b++) acc += D[a][b] * g.Fab[cc][b];
                tmp[a][cc] = acc;
            }
        for (int a = 0; a < 4; a++)
            for (int d = 0; d < 4; d++) {
                Dual acc;
                for (int cc = 0; cc < 4; cc++) acc += D[d][cc] * tmp[a][cc];
                t2[a][d] = -(acc * ipn);
            }
    }
    // t_abc
    Dual q[4], pdn[4];   // q_c = dn[c][b] p^b, pdn_b = p^c dn[c][b]
    for (int a = 0; a < 4; a++) {
        q[a] = Dual(); pdn[a] = Dual();
        for (int b = 0; b < 4; b++) { q[a] += g.Fab[a][b] * p[b]; pdn[a] += p[b] * g.Fab[b][a]; }
    }
    Dual pdnD[4], qD[4];
    for (int a = 0; a < 4; a++) {
        pdnD[a] = Dual(); qD[a] = Dual();
        for (int b = 0; b < 4; b++) { pdnD[a] += pdn[b] * D[a][b]; qD[a] += q[b] * D[a][b]; }
    }
    Dual t3[4][4][4];
    {
        Dual X[4][4][4], Y[4][4][4];   // X[c][e][a] = sum_b ddn[c][e][b] D[a][b]; Y[c][f][a] = sum_e X[c][e][a] D[f][e]
        for (int cc = 0; cc < 4; cc++)
            for (int e = 0; e < 4; e++)
                for (int a = 0; a < 4; a++) {
                    Dual acc;
                    for (int b = 0; b < 4; b++) acc += g.Fabc[cc][e][b] * D[a][b];
                    X[cc][e][a] = acc;
                }
        for (int cc = 0; cc < 4; cc++)
            for (int f = 0; f < 4; f++)
                for (int a = 0; a < 4; a++) {
                    Dual acc;
                    for (int e = 0; e < 4; e++) acc += X[cc][e][a] * D[f][e];
                    Y[cc][f][a] = acc;
                }
        for (int a = 0; a < 4; a++)
            for (int d = 0; d < 4; d++)
                for (int f = 0; f < 4; f++) {
                    Dual acc;
                    for (int cc = 0; cc < 4; cc++) acc += Y[cc][f][a] * D[d][cc];
                    acc += pdnD[a] * t2[d][f] + qD[d] * t2[a][f] + qD[f] * t2[a][d];
                    t3[a][d][f] = -(acc * ipn);
                }
    }
    // Gamma and its derivatives
    Dual QQ[4][4];
    for (int b = 0; b < 4; b++)
        for (int cc = 0; cc < 4; cc++) {
            Dual acc;
            for (int e = 0; e < 4; e++) acc += p[e] * g.Fabc[b][cc][e];
            QQ[b][cc] = acc;
        }
    Dual ipn2 = ipn * ipn, ipn3 = ipn2 * ipn;
    Dual Gam[4], G1[4][4], G2[4][4][4];
    for (int nu = 0; nu < 4; nu++) Gam[nu] = s * n[nu] * ipn;
    for (int b = 0; b < 4; b++)
        for (int nu = 0; nu < 4; nu++) G1[b][nu] = s * (g.Fab[b][nu] * ipn - n[nu] * q[b] * ipn2);
    for (int b = 0; b < 4; b++)
        for (int cc = 0; cc < 4; cc++)
            for (int nu = 0; nu < 4; nu++)
                G2[b][cc][nu] = s * (g.Fabc[b][cc][nu] * ipn
                                     - (g.Fab[b][nu] * q[cc] + g.Fab[cc][nu] * q[b]) * ipn2
                                     - n[nu] * QQ[b][cc] * ipn2
                                     + 2.0 * n[nu] * q[b] * q[cc] * ipn3);
    // raised t's, projected derivatives of Gamma
    Dual t2u[4][4], t3u[4][4][4];
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++) {
            t2u[a][b] = Gm[a] * Gm[b] * t2[a][b];
            for (int cc = 0; cc < 4; cc++) t3u[a][b][cc] = Gm[a] * Gm[b] * Gm[cc] * t3[a][b][cc];
        }
    Dual DG1[4][4], DDG2[4][4][4], pG1[4];
    for (int a = 0; a < 4; a++)
        for (int nu = 0; nu < 4; nu++) {
            Dual acc;
            for (int e = 0; e < 4; e++) acc += Dup[a][e] * G1[e][nu];
            DG1[a][nu] = acc;
        }
    {
        Dual H[4][4][4];  // H[a][f][nu] = sum_e Dup[a][e] G2[e][f][nu]
        for (int a = 0; a < 4; a++)
            for (int f = 0; f < 4; f++)
                for (int nu = 0; nu < 4; nu++) {
                    Dual acc;
                    for (int e = 0; e < 4; e++) acc += Dup[a][e] * G2[e][f][nu];
                    H[a][f][nu] = acc;
                }
        for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
                for (int nu = 0; nu < 4; nu++) {
                    Dual acc;
                    for (int f = 0; f < 4; f++) acc += Dup[b][f] * H[a][f][nu];
                    DDG2[a][b][nu] = acc;
                }
    }
    for (int nu = 0; nu < 4; nu++) {
        pG1[nu] = Dual();
        for (int e = 0; e < 4; e++) pG1[nu] += p[e] * G1[e][nu];
    }
    // contractions with the gradients:  w_s = p^alpha d_s beta_alpha  (d zeta = 0)
    Dual w[4], Dw[4], pdb[4], Ddb[4][4], wp, pGam;
    for (int sg = 0; sg < 4; sg++) {
        w[sg] = Dual();
        for (int al = 0; al < 4; al++) w[sg] += p[al] * c.dbeta[sg][al];
    }
    for (int r = 0; r < 4; r++) {
        pdb[r] = Dual();
        for (int sg = 0; sg < 4; sg++) pdb[r] += p[sg] * c.dbeta[sg][r];
    }
    for (int a = 0; a < 4; a++) {
        Dw[a] = Dual();
        for (int sg = 0; sg < 4; sg++) Dw[a] += Dup[a][sg] * w[sg];
        for (int r = 0; r < 4; r++) {
            Dual acc;
            for (int sg = 0; sg < 4; sg++) acc += Dup[a][sg] * c.dbeta[sg][r];
            Ddb[a][r] = acc;
        }
    }
    wp = Dual(); pGam = Dual();
    for (int a = 0; a < 4; a++) { wp += p[a] * w[a]; pGam += p[a] * Gam[a]; }

    // ---- N = 1  (eq. 46)
    Dual X1[4][4];
    for (int r = 0; r < 4; r++)
        for (int t = 0; t < 4; t++) X1[r][t] = (Gam[r] * Dw[t] + pGam * Ddb[t][r]) * Gm[t];
    eps_contract(pl, X1, R.A1);
    for (int m = 0; m < 4; m++) R.A1[m] = -(a0 * R.A1[m]);

    // ---- N = 3
    Dual pDG1[4], pGY2[4][4], pDDG2[4][4], GY2[4][4][4];
    for (int a = 0; a < 4; a++) {
        pDG1[a] = Dual();
        for (int nu = 0; nu < 4; nu++) pDG1[a] += p[nu] * DG1[a][nu];
        for (int b = 0; b < 4; b++) {
            pGY2[a][b] = Dual(); pDDG2[a][b] = Dual();
            for (int nu = 0; nu < 4; nu++) {
                GY2[a][b][nu] = DDG2[a][b][nu] + pG1[nu] * t2u[a][b];
                pGY2[a][b] += p[nu] * GY2[a][b][nu];
                pDDG2[a][b] += p[nu] * DDG2[a][b][nu];
            }
        }
    }
    // U[r][a][b][c] = V_r^nu_sigma E^{abc sigma}_nu  (complete and literal)
    Dual XB[4][4], XC[4][4], XBl[4][4], XCl[4][4];
    for (int r = 0; r < 4; r++)
        for (int t = 0; t < 4; t++) { XB[r][t] = XC[r][t] = XBl[r][t] = XCl[r][t] = Dual(); }
    for (int r = 0; r < 4; r++) {
        Dual U[4][4][4], Ul[4][4][4];
        for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
                for (int cc = 0; cc < 4; cc++) {
                    U[a][b][cc] = wp * (Gam[r] * t3u[a][b][cc] + 3.0 * DG1[a][r] * t2u[b][cc])
                                  + 3.0 * GY2[a][b][r] * Dw[cc]
                                  + pdb[r] * (pGam * t3u[a][b][cc] + 3.0 * pDG1[a] * t2u[b][cc])
                                  + 3.0 * pGY2[a][b] * Ddb[cc][r];
                    Ul[a][b][cc] = 3.0 * DDG2[a][b][r] * Dw[cc] + 3.0 * pDDG2[a][b] * Ddb[cc][r];
                }
        // T^(3) contractions:  B_{abc tau} = beta_a beta_b g_{c tau} + perms,  C = g_ab g_{c tau} + perms
        for (int t = 0; t < 4; t++) {
            Dual sb, sbl, sc, scl;
            for (int a = 0; a < 4; a++)
                for (int b = 0; b < 4; b++) {
                    double bb = beta[a] * beta[b];
                    sb += bb * (U[a][b][t] + U[a][t][b] + U[t][a][b]);
                    sbl += bb * (Ul[a][b][t] + Ul[a][t][b] + Ul[t][a][b]);
                }
            for (int a = 0; a < 4; a++) {
                sc += Gm[a] * (U[a][a][t] + U[a][t][a] + U[t][a][a]);
                scl += Gm[a] * (Ul[a][a][t] + Ul[a][t][a] + Ul[t][a][a]);
            }
            XB[r][t] = Gm[t] * sb; XBl[r][t] = Gm[t] * sbl;
            XC[r][t] = Gm[t] * sc; XCl[r][t] = Gm[t] * scl;
        }
    }
    Dual eB[4], eC[4], eBl[4], eCl[4];
    eps_contract(pl, XB, eB); eps_contract(pl, XC, eC);
    eps_contract(pl, XBl, eBl); eps_contract(pl, XCl, eCl);
    for (int m = 0; m < 4; m++) {
        R.A3d[m] = (2.0 / 6.0) * a2 * eB[m];
        R.A3dp[m] = (0.5 / 6.0) * a0 * eC[m];
        R.A3d_lit[m] = (2.0 / 6.0) * a2 * eBl[m];
        R.A3dp_lit[m] = (0.5 / 6.0) * a0 * eCl[m];
    }
}

}  // namespace

// ------------------------------------------------------------------ geometry file
void read_geometry(const string &filename, vector<element> &hypersurface, vector<geo_cell> &geometry) {
    ifstream f(filename, ios::binary);
    if (!f.is_open()) {
        cout << "Failed to open " << filename << endl;
        exit(1);
    }
    char head[256];
    f.read(head, 256);
    if (!f || strncmp(head, "FOGEO001", 8) != 0) {
        cout << filename << " is not a geometry file written by write_geometry.py" << endl;
        exit(1);
    }
    int64_t hdr[3];
    memcpy(hdr, head + 8, 3 * sizeof(int64_t));
    int64_t ncells = hdr[0], nrec = hdr[1];
    if (nrec != 42) {
        cout << "unexpected record length " << nrec << " in " << filename << endl;
        exit(1);
    }
    if ((size_t)ncells != hypersurface.size()) {
        cout << "geometry file has " << ncells << " cells, the hypersurface " << hypersurface.size()
             << ": they do not belong together" << endl;
        exit(1);
    }
    cout << "Reading geometry from " << filename << " (" << string(head + 32, strnlen(head + 32, 224)) << ")" << endl;
    geometry.resize(ncells);
    vector<double> rec(42);
    for (int64_t i = 0; i < ncells; i++) {
        f.read(reinterpret_cast<char *>(rec.data()), 42 * sizeof(double));
        if (!f) {
            cout << "geometry file truncated at cell " << i << endl;
            exit(1);
        }
        element &c = hypersurface[i];
        const double pos[4] = {c.tau, c.x, c.y, c.eta};
        for (int k = 0; k < 4; k++)
            if (fabs(rec[k] - pos[k]) > 1e-9 * max(1.0, fabs(pos[k]))) {
                cout << "cell " << i << ": position in the geometry file differs from the hypersurface file" << endl;
                exit(1);
            }
        for (int k = 0; k < 4; k++) c.dsigma[k] = rec[4 + k];      // dSigma aligned with the analytic normal
        geo_cell &g = geometry[i];
        for (int k = 0; k < 4; k++) g.Fa[k] = rec[8 + k];
        for (int j = 0; j < 10; j++) {
            int a = PAIRS[j][0], b = PAIRS[j][1];
            g.Fab[a][b] = g.Fab[b][a] = rec[12 + j];
        }
        int j = 0;
        for (int a = 0; a < 4; a++)
            for (int b = a; b < 4; b++)
                for (int cc = b; cc < 4; cc++, j++) {
                    double v = rec[22 + j];
                    g.Fabc[a][b][cc] = g.Fabc[a][cc][b] = g.Fabc[b][a][cc] = v;
                    g.Fabc[b][cc][a] = g.Fabc[cc][a][b] = g.Fabc[cc][b][a] = v;
                }
    }
    cout << "Geometry read: " << ncells << " cells" << endl;
}

// ------------------------------------------------------------------ polarization
void curvature_polarization(double pT, double phi, double y_rap, pdg_particle particle,
                            vector<element> &hypersurface, vector<geo_cell> &geometry,
                            double cut, bool theta, ofstream &fileout) {
    const double mass = particle.get_mass();
    const int B = particle.get_b(), Q = particle.get_q(), S = particle.get_s();
    const double mT = sqrt(mass * mass + pT * pT);
    const Dual p[4] = {Dual(mT * cosh(y_rap), 1.0), Dual(pT * cos(phi)), Dual(pT * sin(phi)), Dual(mT * sinh(y_rap))};

    double Np = 0;
    double S1[4] = {0, 0, 0, 0}, S3d[4] = {0, 0, 0, 0}, S3dp[4] = {0, 0, 0, 0};
    double S3dl[4] = {0, 0, 0, 0}, S3dpl[4] = {0, 0, 0, 0};
    const long ncell = (long)hypersurface.size();
#ifdef OPEN_MP
#pragma omp parallel for schedule(dynamic, 2048) reduction(+ : Np, S1, S3d, S3dp, S3dl, S3dpl)
#endif
    for (long i = 0; i < ncell; i++) {
        const element &c = hypersurface[i];
        const double mutot = c.mub * B + c.muq * Q + c.mus * S;
        CellResult R;
        cell_terms(p, c, geometry[i], mass, mutot, R);
        if (theta && R.pdS.v < 0) continue;                              // as --impr-theta
        Np += R.pdS.v * R.nF.v;
        for (int m = 0; m < 4; m++) S1[m] += R.pdS.v * R.A1[m].v;     // eq. (46): all cells
        if (R.grazing < cut) continue;                                   // cut acts on the N = 3 term only
        for (int m = 0; m < 4; m++) {
            S3d[m] += R.pdS.v * R.A3d[m].v;
            S3dl[m] += R.pdS.v * R.A3d_lit[m].v;
            // delta' term: - d/dp0 [ (dSigma.p) A_delta' / (2 p0) ]
            Dual gq = R.pdS * R.A3dp[m] / (2.0 * p[0]);
            Dual gl = R.pdS * R.A3dp_lit[m] / (2.0 * p[0]);
            S3dp[m] += -gq.d;
            S3dpl[m] += -gl.d;
        }
    }
    const double n1 = hbarC / (8.0 * mass * Np);
    const double n3 = hbarC * hbarC * hbarC / (8.0 * mass * Np);
    fileout << setprecision(10) << "   " << pT << "   " << phi << "   " << y_rap << "   " << Np;
    for (int m = 0; m < 4; m++) fileout << "   " << S1[m] * n1;
    for (int m = 0; m < 4; m++) fileout << "   " << (S3d[m] + S3dp[m]) * n3;
    for (int m = 0; m < 4; m++) fileout << "   " << (S3dl[m] + S3dpl[m]) * n3;
    fileout << endl;
}
