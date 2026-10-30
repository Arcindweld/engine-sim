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
    m_invDx = 0.0;
    m_invG1 = 1.0 / (m_gamma - 1);
    m_maxSignalSpeed = 0.0;
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
    m_invA.assign(n, 0.0);
    m_invD.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        m_invA[i] = 1.0 / m_A[i];
        m_invD[i] = 1.0 / m_D[i];
    }
    m_invDx = 1.0 / m_dx;
    m_invG1 = 1.0 / (m_gamma - 1);
    m_invR.assign(n, 0.0);

    m_q1.assign(n, 0.0);
    m_q2.assign(n, 0.0);
    m_q3.assign(n, 0.0);
    m_r.assign(n, 0.0);
    m_u.assign(n, 0.0);
    m_p.assign(n, 0.0);
    for (int k = 0; k < ScalarCount; ++k) {
        m_qs[k].assign(n, 0.0);
        m_x[k].assign(n, 0.0);
    }
    for (int k = 0; k < FluxCount; ++k) m_f[k].assign(n + 1, 0.0);
    m_bL.assign(n, State());
    m_bR.assign(n, State());
    m_maxSignalSpeed = 0.0;

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
    double smax = 0.0;
    for (int i = 0; i < m_n; ++i) {
        const double rr = m_q1[i] / m_A[i];
        const double uu = m_q2[i] / m_q1[i];
        const double E = m_q3[i] / m_A[i];
        m_r[i] = rr;
        m_invR[i] = 1.0 / rr;
        m_u[i] = uu;
        m_p[i] = std::max(1.0, g1 * (E - 0.5 * rr * uu * uu));
        for (int k = 0; k < ScalarCount; ++k) {
            m_x[k][i] = std::min(1.0, std::max(0.0, m_qs[k][i] / m_q1[i]));
        }
        const double v = std::fabs(uu) + std::sqrt(m_gamma * m_p[i] / rr);
        if (v > smax) smax = v;
    }
    m_maxSignalSpeed = smax;
}

Pipe1D::State Pipe1D::state(int i) const {
    State s;
    s.rho = m_r[i];
    s.u = m_u[i];
    s.p = m_p[i];
    for (int k = 0; k < ScalarCount; ++k) s.x[k] = m_x[k][i];
    return s;
}

// HLLC flux per unit area (Toro, 3rd ed., sec. 10.4-10.5), with the passive
// scalars carried by the mass flux of the side holding the contact wave.
void Pipe1D::hllc(const State &L, const State &R, double *F) const {
    const double invG1 = m_invG1;
    const double rL = L.rho, uL = L.u, pL = L.p;
    const double rR = R.rho, uR = R.u, pR = R.p;
    const double invRL = 1.0 / rL, invRR = 1.0 / rR;
    const double cL = std::sqrt(m_gamma * pL * invRL), cR = std::sqrt(m_gamma * pR * invRR);
    const double EL = pL * invG1 + 0.5 * rL * uL * uL, ER = pR * invG1 + 0.5 * rR * uR * uR;
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
        const double w = SL - uL;
        const double k = rL * w / (SL - Ss);
        F[0] = rL * uL + SL * (k - rL);
        F[1] = rL * uL * uL + pL + SL * (k * Ss - rL * uL);
        F[2] = uL * (EL + pL) + SL * (k * (EL * invRL + (Ss - uL) * (Ss + pL * invRL / w)) - EL);
        for (int s = 0; s < ScalarCount; ++s) F[3 + s] = F[0] * L.x[s];
    }
    else {
        const double w = SR - uR;
        const double k = rR * w / (SR - Ss);
        F[0] = rR * uR + SR * (k - rR);
        F[1] = rR * uR * uR + pR + SR * (k * Ss - rR * uR);
        F[2] = uR * (ER + pR) + SR * (k * (ER * invRR + (Ss - uR) * (Ss + pR * invRR / w)) - ER);
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

// Open end to still gas at (pAtm, TAtm): method-of-characteristics boundary
// (Blair, "Design and Simulation of Four-Stroke Engines", sec. 2.16-2.17).
// The Riemann invariant leaving the pipe through the end cell is kept,
// J = u + 2c/(gamma-1) (right end) or u - 2c/(gamma-1) (left end).
//   outflow: boundary pressure = pAtm, state reached isentropically from the
//            end cell; a supersonic exit keeps the cell state.
//   inflow:  the gas comes from rest at (pAtm, TAtm) and expands isentropically
//            to the boundary pressure, cp T_b + u_b^2/2 = cp TAtm,
//            p_b = pAtm (T_b/TAtm)^(gamma/(gamma-1)), while the invariant of
//            the (possibly much hotter) pipe gas fixes u_b for that p_b: the
//            contact between the two gases sits at the boundary, pressure and
//            velocity are continuous across it, the sound speed is not.  For a
//            small suction this gives the acoustic inflow u ~ dp/(rho c), which
//            keeps a pipe at rest at rest (a Bernoulli ghost velocity ~ sqrt(dp)
//            is unstable to round-off).
void Pipe1D::applyOpenEnd(int side, double pAtm, double TAtm, const double *scalarsAtm) {
    const int i = endCell(side), j = endFace(side);
    const double g = m_gamma, g1 = g - 1, k2 = 2 / g1;
    const double r = m_r[i], u = m_u[i], p = m_p[i];
    const double c = std::sqrt(g * p / r);
    const double af = m_Af[j];

    // outward-facing velocity and the invariant carried out of the pipe
    const double un = u * side;
    const double J = un + k2 * c;

    double rb, ub, pb;               // boundary state, ub in the outward direction
    bool inflow = false;
    if (un >= c && p >= pAtm) {
        // supersonic exit above the outside pressure: nothing comes back
        rb = r; ub = un; pb = p;
    }
    else {
        const double cb = c * std::pow(pAtm / p, g1 / (2 * g));
        const double ubOut = J - k2 * cb;
        if (ubOut >= 0) {
            pb = pAtm;
            rb = r * std::pow(pAtm / p, 1 / g);
            ub = ubOut;
        }
        else {
            inflow = true;
            // solve  v - J + k2 * c * (pb(v) / p)^ex = 0  for v < 0, with the
            // cold side's Bernoulli pressure pb(v) = pAtm (1 - v^2/(2 cp TAtm))^(g/g1)
            const double ex = g1 / (2 * g);
            const double ca2 = g * (m_R * TAtm);
            const double vMax = std::sqrt(2 * ca2 / (g + 1));   // sonic inflow from stagnation
            const double invCpTa = 1.0 / (m_cp * TAtm);
            double v = ubOut;
            if (-v > vMax) v = -vMax;
            for (int it = 0; it < 30; ++it) {
                const double tb = std::max(1e-6, 1 - 0.5 * v * v * invCpTa);   // T_b / TAtm
                const double pbv = pAtm * std::pow(tb, g / g1);
                const double ratio = std::pow(pbv / p, ex);
                const double f = v - J + k2 * c * ratio;
                const double dpb = pAtm * (g / g1) * std::pow(tb, g / g1 - 1) * (-v * invCpTa);
                const double df = 1 + k2 * c * ex * ratio / pbv * dpb;
                const double step = f / df;
                v -= step;
                if (v > 0) v = 0;
                if (-v > vMax) v = -vMax;
                if (std::fabs(step) < 1e-9 * std::max(1.0, std::fabs(v))) break;
            }
            const double Tb = TAtm - 0.5 * v * v / m_cp;
            pb = pAtm * std::pow(Tb / TAtm, g / g1);
            rb = pb / (m_R * Tb);
            ub = v;
        }
    }

    // Euler flux of the boundary state in the +x direction
    const double ux = ub * side;
    const double Eb = pb / g1 + 0.5 * rb * ub * ub;
    m_f[0][j] = rb * ux * af;
    m_f[1][j] = (rb * ux * ux + pb) * af;
    m_f[2][j] = ux * (Eb + pb) * af;
    for (int k = 0; k < ScalarCount; ++k) {
        const double xk = inflow ? scalarsAtm[k] : m_x[k][i];
        m_f[3 + k][j] = rb * ux * xk * af;
    }
}

// One MUSCL-Hancock step (reference: stepPipe).  Slopes are minmod-limited
// differences of the primitive variables; the half-step predictor uses the
// linearised primitive equations; interior face fluxes come from HLLC; the
// end fluxes were set by the owner.  Sources: area change p dA/dx, wall
// friction 2 f rho |u| A / D and Reynolds-analogy heat transfer.
void Pipe1D::advance(double dt) {
    const int n = m_n;
    const double g1 = m_gamma - 1;
    const double gamma = m_gamma;
    const double h = 0.5 * dt / m_dx;

    // raw pointers: lets the compiler keep the array bases in registers
    const double *A = m_A.data(), *Af = m_Af.data(), *invA = m_invA.data(), *invD = m_invD.data();
    double *q1 = m_q1.data(), *q2 = m_q2.data(), *q3 = m_q3.data();
    double *r = m_r.data(), *u = m_u.data(), *p = m_p.data(), *invR = m_invR.data();
    double *qs[ScalarCount], *x[ScalarCount], *f[FluxCount];
    for (int k = 0; k < ScalarCount; ++k) { qs[k] = m_qs[k].data(); x[k] = m_x[k].data(); }
    for (int k = 0; k < FluxCount; ++k) f[k] = m_f[k].data();
    State *bL = m_bL.data(), *bR = m_bR.data();

    // Limited slopes and the half-step predictor, once per cell: the evolved
    // left (bL) and right (bR) face states of each cell.
    for (int i = 0; i < n; ++i) {
        double dr = 0.0, du = 0.0, dp = 0.0, dxs[ScalarCount] = { 0.0, 0.0 };
        if (i > 0 && i < n - 1) {
            dr = minmod(r[i] - r[i - 1], r[i + 1] - r[i]);
            du = minmod(u[i] - u[i - 1], u[i + 1] - u[i]);
            dp = minmod(p[i] - p[i - 1], p[i + 1] - p[i]);
            for (int k = 0; k < ScalarCount; ++k) {
                dxs[k] = minmod(x[k][i] - x[k][i - 1], x[k][i + 1] - x[k][i]);
            }
        }
        const double ri = r[i], ui = u[i], pi = p[i];
        const double er = -h * (ui * dr + ri * du);
        const double eu = -h * (ui * du + dp * invR[i]);
        const double ep = -h * (gamma * pi * du + ui * dp);
        State &L = bL[i];
        State &R = bR[i];
        L.rho = ri - 0.5 * dr + er; R.rho = ri + 0.5 * dr + er;
        L.u = ui - 0.5 * du + eu;   R.u = ui + 0.5 * du + eu;
        L.p = pi - 0.5 * dp + ep;   R.p = pi + 0.5 * dp + ep;
        for (int k = 0; k < ScalarCount; ++k) {
            const double ex = -h * ui * dxs[k];
            L.x[k] = x[k][i] - 0.5 * dxs[k] + ex;
            R.x[k] = x[k][i] + 0.5 * dxs[k] + ex;
        }
        if (!(L.rho > 0 && L.p > 0)) L = state(i);
        if (!(R.rho > 0 && R.p > 0)) R = state(i);
    }

    double F[FluxCount];
    for (int j = 1; j < n; ++j) {
        hllc(bR[j - 1], bL[j], F);
        const double af = Af[j];
        for (int k = 0; k < FluxCount; ++k) f[k][j] = F[k] * af;
    }

    // Conservative update with the sources, then the new primitive state and
    // the largest signal speed for the owner's CFL estimate.
    const double kdt = dt * m_invDx;
    const double invRgas = 1.0 / m_R;
    const double friction2 = 2 * m_friction;
    const double cpTw = m_cp, Tw = m_Tw;
    double smax = 0.0;
    for (int i = 0; i < n; ++i) {
        const double ui = u[i], ri = r[i], pi = p[i], invRi = invR[i];
        const double Ti = pi * invRgas * invRi;
        const double fr = friction2 * ri * std::fabs(ui) * A[i] * invD[i];

        double m1 = q1[i] - kdt * (f[0][i + 1] - f[0][i]);
        const double m2 = q2[i] - kdt * (f[1][i + 1] - f[1][i])
            + dt * (pi * (Af[i + 1] - Af[i]) * m_invDx - fr * ui);
        const double m3 = q3[i] - kdt * (f[2][i + 1] - f[2][i]) + dt * fr * cpTw * (Tw - Ti);
        const double floorMass = MinDensityFraction * A[i];
        if (m1 < floorMass) m1 = floorMass;
        q1[i] = m1;
        q2[i] = m2;
        q3[i] = m3;
        const double invQ1 = 1.0 / m1;
        for (int k = 0; k < ScalarCount; ++k) {
            double v = qs[k][i] - kdt * (f[3 + k][i + 1] - f[3 + k][i]);
            if (v < 0) v = 0;
            else if (v > m1) v = m1;
            qs[k][i] = v;
            x[k][i] = v * invQ1;
        }

        const double rr = m1 * invA[i];
        const double uu = m2 * invQ1;
        const double E = m3 * invA[i];
        const double pp = std::max(1.0, g1 * (E - 0.5 * rr * uu * uu));
        const double ir = A[i] * invQ1;
        r[i] = rr;
        invR[i] = ir;
        u[i] = uu;
        p[i] = pp;
        const double v = std::fabs(uu) + std::sqrt(gamma * pp * ir);
        if (v > smax) smax = v;
    }
    m_maxSignalSpeed = smax;
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
    m_invR[i] = 1.0 / rho;
    m_u[i] = u;
    m_p[i] = p;
    for (int k = 0; k < ScalarCount; ++k) m_x[k][i] = scalars ? scalars[k] : 0.0;
}

void Pipe1D::addToCell(int i, double dMass, double dMomentum, double dEnergy, const double *dScalarMass) {
    const double inv = 1.0 / m_dx;
    m_q1[i] += dMass * inv;
    m_q2[i] += dMomentum * inv;
    m_q3[i] += dEnergy * inv;
    if (m_q1[i] < MinDensityFraction * m_A[i]) m_q1[i] = MinDensityFraction * m_A[i];
    for (int k = 0; k < ScalarCount; ++k) {
        m_qs[k][i] += dScalarMass[k] * inv;
        if (m_qs[k][i] < 0) m_qs[k][i] = 0;
        else if (m_qs[k][i] > m_q1[i]) m_qs[k][i] = m_q1[i];
    }
    const double rr = m_q1[i] * m_invA[i];
    const double uu = m_q2[i] / m_q1[i];
    const double E = m_q3[i] * m_invA[i];
    m_r[i] = rr;
    m_invR[i] = 1.0 / rr;
    m_u[i] = uu;
    m_p[i] = std::max(1.0, (m_gamma - 1) * (E - 0.5 * rr * uu * uu));
    for (int k = 0; k < ScalarCount; ++k) m_x[k][i] = m_qs[k][i] / m_q1[i];
    const double v = std::fabs(uu) + std::sqrt(m_gamma * m_p[i] * m_invR[i]);
    if (v > m_maxSignalSpeed) m_maxSignalSpeed = v;
}
