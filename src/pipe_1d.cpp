// One-dimensional compressible duct solver for exhaust pipes.
//
// Copyright (c) 2026 Arcindweld. MIT License (see LICENSE).
//
// See pipe_1d.h for the model and the reference it was ported from.

#include "../include/pipe_1d.h"

#include <algorithm>
#include <cmath>

namespace {
    constexpr double Pi = 3.14159265358979323846;
    constexpr double MinDensityFraction = 1e-6;   // floor on rho A, as in the reference
}

Pipe1D::Pipe1D() {
    m_n = 0;
    m_dx = 0.0;
    m_gamma = 1.4;
    m_R = 287.0;
    m_cv = m_R / (m_gamma - 1);
    m_cp = m_gamma * m_cv;
    m_Tw = 700.0;
    m_friction = 0.005;
    m_avgMassFlow[0] = m_avgMassFlow[1] = 0.0;
}

Pipe1D::~Pipe1D() {
    /* void */
}

void Pipe1D::initialize(const Parameters &params) {
    m_gamma = params.heatCapacityRatio;
    m_R = params.specificGasConstant;
    m_cv = m_R / (m_gamma - 1);
    m_cp = m_gamma * m_cv;
    m_Tw = params.wallTemperature;
    m_friction = params.frictionFactor;

    double L = 0.0;
    for (const Section &s : params.sections) L += s.length;

    int n = static_cast<int>(std::lround(L / params.cellSize));
    n = std::max(n, params.minCells);
    n = std::min(n, std::max(params.maxCells, params.minCells));
    m_n = n;
    m_dx = L / n;

    // area(x): linear taper of the diameter inside each section (cone)
    auto area = [&](double x) {
        double a = 0.0;
        const size_t count = params.sections.size();
        for (size_t s = 0; s < count; ++s) {
            const Section &sec = params.sections[s];
            if (x <= a + sec.length || s + 1 == count) {
                const double f = (sec.length > 0)
                    ? std::min(1.0, std::max(0.0, (x - a) / sec.length))
                    : 1.0;
                const double d = sec.inletDiameter + (sec.outletDiameter - sec.inletDiameter) * f;
                return Pi * d * d / 4;
            }
            a += sec.length;
        }
        return 0.0;
    };

    m_A.assign(n, 0.0);
    m_Af.assign(n + 1, 0.0);
    m_D.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        m_A[i] = area((i + 0.5) * m_dx);
        m_D[i] = std::sqrt(4 * m_A[i] / Pi);
    }
    m_Af[0] = m_A[0];
    m_Af[n] = m_A[n - 1];
    for (int j = 1; j < n; ++j) m_Af[j] = 0.5 * (m_A[j - 1] + m_A[j]);

    m_q1.assign(n, 0.0);
    m_q2.assign(n, 0.0);
    m_q3.assign(n, 0.0);
    m_r.assign(n, 0.0);
    m_u.assign(n, 0.0);
    m_p.assign(n, 0.0);
    m_dr.assign(n, 0.0);
    m_du.assign(n, 0.0);
    m_dp.assign(n, 0.0);
    for (int k = 0; k < ScalarCount; ++k) {
        m_qs[k].assign(n, 0.0);
        m_x[k].assign(n, 0.0);
        m_dxs[k].assign(n, 0.0);
    }
    for (int k = 0; k < FluxCount; ++k) m_f[k].assign(n + 1, 0.0);

    const double rho0 = params.initialPressure / (m_R * params.initialTemperature);
    for (int i = 0; i < n; ++i) {
        m_q1[i] = rho0 * m_A[i];
        m_q2[i] = 0.0;
        m_q3[i] = params.initialPressure / (m_gamma - 1) * m_A[i];
        for (int k = 0; k < ScalarCount; ++k) m_qs[k][i] = rho0 * params.initialScalars[k] * m_A[i];
    }

    m_avgMassFlow[0] = m_avgMassFlow[1] = 0.0;
    computePrimitives();
}

void Pipe1D::computePrimitives() {
    const double g1 = m_gamma - 1;
    for (int i = 0; i < m_n; ++i) {
        const double rr = m_q1[i] / m_A[i];
        const double uu = m_q2[i] / m_q1[i];
        const double E = m_q3[i] / m_A[i];
        m_r[i] = rr;
        m_u[i] = uu;
        m_p[i] = std::max(1.0, g1 * (E - 0.5 * rr * uu * uu));
        for (int k = 0; k < ScalarCount; ++k) {
            m_x[k][i] = std::min(1.0, std::max(0.0, m_qs[k][i] / m_q1[i]));
        }
    }
}

Pipe1D::State Pipe1D::state(int i) const {
    State s;
    s.rho = m_r[i];
    s.u = m_u[i];
    s.p = m_p[i];
    for (int k = 0; k < ScalarCount; ++k) s.x[k] = m_x[k][i];
    return s;
}

double Pipe1D::maxSignalSpeed() const {
    double s = 0.0;
    for (int i = 0; i < m_n; ++i) {
        const double v = std::fabs(m_u[i]) + std::sqrt(m_gamma * m_p[i] / m_r[i]);
        if (v > s) s = v;
    }
    return s;
}

// HLLC flux per unit area (Toro, 3rd ed., sec. 10.4-10.5), with the passive
// scalars carried by the mass flux of the side holding the contact wave.
void Pipe1D::hllc(const State &L, const State &R, double *F) const {
    const double g1 = m_gamma - 1;
    const double rL = L.rho, uL = L.u, pL = L.p;
    const double rR = R.rho, uR = R.u, pR = R.p;
    const double cL = std::sqrt(m_gamma * pL / rL), cR = std::sqrt(m_gamma * pR / rR);
    const double EL = pL / g1 + 0.5 * rL * uL * uL, ER = pR / g1 + 0.5 * rR * uR * uR;
    const double SL = std::min(uL - cL, uR - cR), SR = std::max(uL + cL, uR + cR);

    if (SL >= 0) {
        F[0] = rL * uL;
        F[1] = rL * uL * uL + pL;
        F[2] = uL * (EL + pL);
        for (int k = 0; k < ScalarCount; ++k) F[3 + k] = F[0] * L.x[k];
        return;
    }
    if (SR <= 0) {
        F[0] = rR * uR;
        F[1] = rR * uR * uR + pR;
        F[2] = uR * (ER + pR);
        for (int k = 0; k < ScalarCount; ++k) F[3 + k] = F[0] * R.x[k];
        return;
    }

    const double Ss =
        (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR))
        / (rL * (SL - uL) - rR * (SR - uR));
    if (Ss >= 0) {
        const double k = rL * (SL - uL) / (SL - Ss);
        F[0] = rL * uL + SL * (k - rL);
        F[1] = rL * uL * uL + pL + SL * (k * Ss - rL * uL);
        F[2] = uL * (EL + pL) + SL * (k * (EL / rL + (Ss - uL) * (Ss + pL / (rL * (SL - uL)))) - EL);
        for (int s = 0; s < ScalarCount; ++s) F[3 + s] = F[0] * L.x[s];
    }
    else {
        const double k = rR * (SR - uR) / (SR - Ss);
        F[0] = rR * uR + SR * (k - rR);
        F[1] = rR * uR * uR + pR + SR * (k * Ss - rR * uR);
        F[2] = uR * (ER + pR) + SR * (k * (ER / rR + (Ss - uR) * (Ss + pR / (rR * (SR - uR)))) - ER);
        for (int s = 0; s < ScalarCount; ++s) F[3 + s] = F[0] * R.x[s];
    }
}

void Pipe1D::setEndFlux(int side, const double *flux) {
    const int j = endFace(side);
    for (int k = 0; k < FluxCount; ++k) m_f[k][j] = flux[k];
}

void Pipe1D::getEndFlux(int side, double *flux) const {
    const int j = endFace(side);
    for (int k = 0; k < FluxCount; ++k) flux[k] = m_f[k][j];
}

void Pipe1D::applyWallEnd(int side) {
    const int i = endCell(side), j = endFace(side);
    m_f[0][j] = 0.0;
    m_f[1][j] = m_p[i] * m_Af[j];
    m_f[2][j] = 0.0;
    for (int k = 0; k < ScalarCount; ++k) m_f[3 + k][j] = 0.0;
}

// Ghost state outside an open end (reference: openGhost).  dir = side, i.e.
// +1 when the end faces +x.  Outflow: the outside sets the pressure unless the
// exit is supersonic and above the outside pressure.  Inflow: still gas at
// (pAtm, TAtm) expanding isentropically to the end cell's pressure.
Pipe1D::State Pipe1D::openEndGhost(int side, double pAtm, double TAtm, const double *scalarsAtm) const {
    const int i = endCell(side);
    const double r = m_r[i], u = m_u[i], p = m_p[i];
    const double g1 = m_gamma - 1;
    State g;
    const double un = u * side;
    if (un >= 0) {
        if (un >= std::sqrt(m_gamma * p / r) && p >= pAtm) {
            g.rho = r; g.u = u; g.p = p;
            for (int k = 0; k < ScalarCount; ++k) g.x[k] = m_x[k][i];
            return g;
        }
        g.p = pAtm;
        g.rho = r * std::pow(pAtm / p, 1 / m_gamma);
        g.u = u;
        for (int k = 0; k < ScalarCount; ++k) g.x[k] = m_x[k][i];
        return g;
    }
    const double pg = std::min(p, pAtm);
    const double Tg = TAtm * std::pow(pg / pAtm, g1 / m_gamma);
    const double v = std::sqrt(std::max(0.0, 2 * m_cp * (TAtm - Tg)));
    g.p = pg;
    g.rho = pg / (m_R * Tg);
    g.u = -side * v;
    for (int k = 0; k < ScalarCount; ++k) g.x[k] = scalarsAtm[k];
    return g;
}

void Pipe1D::applyOpenEnd(int side, double pAtm, double TAtm, const double *scalarsAtm) {
    const int i = endCell(side), j = endFace(side);
    const State g = openEndGhost(side, pAtm, TAtm, scalarsAtm);
    const State c = state(i);
    double F[FluxCount];
    if (side > 0) hllc(c, g, F);
    else hllc(g, c, F);
    const double af = m_Af[j];
    for (int k = 0; k < FluxCount; ++k) m_f[k][j] = F[k] * af;
}

// One MUSCL-Hancock step (reference: stepPipe).  Slopes are minmod-limited
// differences of the primitive variables; the half-step predictor uses the
// linearised primitive equations; interior face fluxes come from HLLC; the
// end fluxes were set by the owner.  Sources: area change p dA/dx, wall
// friction 2 f rho |u| A / D and Reynolds-analogy heat transfer.
void Pipe1D::advance(double dt) {
    const int n = m_n;
    const double g1 = m_gamma - 1;

    for (int i = 0; i < n; ++i) {
        if (i == 0 || i == n - 1) {
            m_dr[i] = m_du[i] = m_dp[i] = 0.0;
            for (int k = 0; k < ScalarCount; ++k) m_dxs[k][i] = 0.0;
            continue;
        }
        m_dr[i] = minmod(m_r[i] - m_r[i - 1], m_r[i + 1] - m_r[i]);
        m_du[i] = minmod(m_u[i] - m_u[i - 1], m_u[i + 1] - m_u[i]);
        m_dp[i] = minmod(m_p[i] - m_p[i - 1], m_p[i + 1] - m_p[i]);
        for (int k = 0; k < ScalarCount; ++k) {
            m_dxs[k][i] = minmod(m_x[k][i] - m_x[k][i - 1], m_x[k][i + 1] - m_x[k][i]);
        }
    }

    const double h = 0.5 * dt / m_dx;
    State a, b;
    double F[FluxCount];
    for (int j = 1; j < n; ++j) {
        const int ia = j - 1, ib = j;

        a.rho = m_r[ia] + 0.5 * m_dr[ia];
        a.u = m_u[ia] + 0.5 * m_du[ia];
        a.p = m_p[ia] + 0.5 * m_dp[ia];
        a.rho -= h * (m_u[ia] * m_dr[ia] + m_r[ia] * m_du[ia]);
        a.u -= h * (m_u[ia] * m_du[ia] + m_dp[ia] / m_r[ia]);
        a.p -= h * (m_gamma * m_p[ia] * m_du[ia] + m_u[ia] * m_dp[ia]);
        for (int k = 0; k < ScalarCount; ++k) {
            a.x[k] = m_x[k][ia] + 0.5 * m_dxs[k][ia] - h * m_u[ia] * m_dxs[k][ia];
        }

        b.rho = m_r[ib] - 0.5 * m_dr[ib];
        b.u = m_u[ib] - 0.5 * m_du[ib];
        b.p = m_p[ib] - 0.5 * m_dp[ib];
        b.rho -= h * (m_u[ib] * m_dr[ib] + m_r[ib] * m_du[ib]);
        b.u -= h * (m_u[ib] * m_du[ib] + m_dp[ib] / m_r[ib]);
        b.p -= h * (m_gamma * m_p[ib] * m_du[ib] + m_u[ib] * m_dp[ib]);
        for (int k = 0; k < ScalarCount; ++k) {
            b.x[k] = m_x[k][ib] - 0.5 * m_dxs[k][ib] - h * m_u[ib] * m_dxs[k][ib];
        }

        if (!(a.rho > 0 && a.p > 0)) a = state(ia);
        if (!(b.rho > 0 && b.p > 0)) b = state(ib);

        hllc(a, b, F);
        const double af = m_Af[j];
        for (int k = 0; k < FluxCount; ++k) m_f[k][j] = F[k] * af;
    }

    const double kdt = dt / m_dx;
    for (int i = 0; i < n; ++i) {
        const double ui = m_u[i], ri = m_r[i];
        const double Ti = m_p[i] / (m_R * ri);
        const double fr = 2 * m_friction * ri * std::fabs(ui) * m_A[i] / m_D[i];

        m_q1[i] -= kdt * (m_f[0][i + 1] - m_f[0][i]);
        m_q2[i] += -kdt * (m_f[1][i + 1] - m_f[1][i])
            + dt * (m_p[i] * (m_Af[i + 1] - m_Af[i]) / m_dx - fr * ui);
        m_q3[i] += -kdt * (m_f[2][i + 1] - m_f[2][i]) + dt * fr * m_cp * (m_Tw - Ti);
        for (int k = 0; k < ScalarCount; ++k) {
            m_qs[k][i] -= kdt * (m_f[3 + k][i + 1] - m_f[3 + k][i]);
        }

        if (m_q1[i] < MinDensityFraction * m_A[i]) m_q1[i] = MinDensityFraction * m_A[i];
        for (int k = 0; k < ScalarCount; ++k) {
            if (m_qs[k][i] < 0) m_qs[k][i] = 0;
            else if (m_qs[k][i] > m_q1[i]) m_qs[k][i] = m_q1[i];
        }
    }
}

double Pipe1D::totalMass() const {
    double m = 0.0;
    for (int i = 0; i < m_n; ++i) m += m_q1[i] * m_dx;
    return m;
}

double Pipe1D::totalEnergy() const {
    double e = 0.0;
    for (int i = 0; i < m_n; ++i) e += m_q3[i] * m_dx;
    return e;
}

void Pipe1D::setCellState(int i, double rho, double u, double p, const double *scalars) {
    m_q1[i] = rho * m_A[i];
    m_q2[i] = rho * u * m_A[i];
    m_q3[i] = (p / (m_gamma - 1) + 0.5 * rho * u * u) * m_A[i];
    for (int k = 0; k < ScalarCount; ++k) m_qs[k][i] = rho * (scalars ? scalars[k] : 0.0) * m_A[i];
    m_r[i] = rho;
    m_u[i] = u;
    m_p[i] = p;
    for (int k = 0; k < ScalarCount; ++k) m_x[k][i] = scalars ? scalars[k] : 0.0;
}
