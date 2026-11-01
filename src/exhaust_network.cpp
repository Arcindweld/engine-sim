// Network of 1-D pipes for one exhaust system.
//
// Copyright (c) 2026 Arcindweld. MIT License (see LICENSE).
//
// See exhaust_network.h.  The junction solvers are ports of `junction` and
// `junctionN` from the ExhaustTuner `wsFactory()` reference: characteristic
// (method of characteristics) matching with one pressure at the junction,
// isentropic change of state along each incoming characteristic, mass flows
// that sum to zero, and the inflows mixed before they are sent down the
// branches that receive gas.

#include "../include/exhaust_network.h"

#include "../include/constants.h"

#include <algorithm>
#include <cmath>

namespace {
    constexpr double Pi = 3.14159265358979323846;
    // Same role as the 0.9 clamp in GasSystem::flow for the port volume.
    constexpr double MaxForwardMoleFraction = 0.9;
}

ExhaustNetwork::ExhaustNetwork() {
    m_built = false;
    m_gamma = 1.4;
    m_R = 287.0;
    m_cp = m_gamma / (m_gamma - 1) * m_R;
    m_prc = std::pow(2 / (m_gamma + 1), m_gamma / (m_gamma - 1));
    m_psic = std::sqrt(m_gamma) * std::pow(2 / (m_gamma + 1), (m_gamma + 1) / (2 * (m_gamma - 1)));
    m_atmScalars[0] = m_atmScalars[1] = 0.0;
    m_outletMassFlow = 0.0;
    m_outletMassFlowDerivative = 0.0;
    m_massIn = m_massOut = m_initialMass = 0.0;
    m_lastSubSteps = 0;
    m_stepOutflow = m_stepTime = 0.0;
    m_stepSubSteps = 0;
}

ExhaustNetwork::~ExhaustNetwork() {
    for (Pipe1D *p : m_collector) delete p;
    m_collector.clear();
}

int ExhaustNetwork::addPrimary(const PrimarySpec &spec) {
    Primary prim;
    prim.spec = spec;
    prim.lastMassFlow = 0.0;
    prim.facePressure = 0.0;
    m_primaries.push_back(prim);
    return static_cast<int>(m_primaries.size()) - 1;
}

void ExhaustNetwork::build(const Parameters &params) {
    m_params = params;
    m_gamma = 1.0 + 2.0 / params.degreesOfFreedom;
    m_R = constants::R / params.molarMass;
    m_cp = m_gamma / (m_gamma - 1) * m_R;
    m_prc = std::pow(2 / (m_gamma + 1), m_gamma / (m_gamma - 1));
    m_psic = std::sqrt(m_gamma) * std::pow(2 / (m_gamma + 1), (m_gamma + 1) / (2 * (m_gamma - 1)));

    // The legacy model treats the air outside the exhaust as inert (ExhaustSystem::process)
    m_atmScalars[0] = 0.0;   // fuel
    m_atmScalars[1] = 0.0;   // O2

    Pipe1D::Parameters pp;
    pp.cellSize = params.cellSize;
    pp.maxCells = params.maxCellsPerPipe;
    pp.wallTemperature = params.wallTemperature;
    pp.frictionFactor = params.frictionFactor;
    pp.initialPressure = params.atmosphericPressure;
    pp.initialTemperature = params.initialTemperature;
    pp.heatCapacityRatio = m_gamma;
    pp.specificGasConstant = m_R;
    pp.initialScalars[0] = m_atmScalars[0];
    pp.initialScalars[1] = m_atmScalars[1];

    for (Primary &prim : m_primaries) {
        Pipe1D::Section s;
        // a primary needs some length to carry cells: a zero length (primary_tube_length 0 with no per-cylinder
        // length) gave a zero cell size, a CFL loop that never advanced and NaN mass. Under a millimetre it
        // becomes one cell length (a port stub); real short primaries are left as given.
        s.length = prim.spec.length < 1e-3 ? params.cellSize : prim.spec.length;
        s.inletDiameter = s.outletDiameter = prim.spec.diameter;
        pp.sections.assign(1, s);
        prim.pipe.initialize(pp);
    }

    // Collector: consecutive sections with matching diameters share a pipe
    // (cones stay inside it); a diameter jump starts a new pipe joined by a
    // characteristic step junction (reference: buildLine).
    for (Pipe1D *p : m_collector) delete p;
    m_collector.clear();
    m_collectorJunctions.clear();

    std::vector<std::vector<Pipe1D::Section>> runs;
    for (const Pipe1D::Section &sec : params.collectorSections) {
        if (!(sec.length > 1e-6)) continue;
        if (!runs.empty() && std::fabs(runs.back().back().outletDiameter - sec.inletDiameter) < 1e-6) {
            runs.back().push_back(sec);
        }
        else {
            runs.push_back(std::vector<Pipe1D::Section>(1, sec));
        }
    }
    if (runs.empty()) {
        // No collector given: a short stub of the primary diameter so that the
        // junction and the open end still exist.
        Pipe1D::Section s;
        s.length = params.cellSize * pp.minCells;
        s.inletDiameter = s.outletDiameter = m_primaries.empty() ? 0.05 : m_primaries[0].spec.diameter;
        runs.push_back(std::vector<Pipe1D::Section>(1, s));
    }
    for (const auto &run : runs) {
        Pipe1D *pipe = new Pipe1D;
        pp.sections = run;
        pipe->initialize(pp);
        m_collector.push_back(pipe);
    }
    m_collectorJunctions.assign(m_collector.size(), Junction());

    m_junctionEnds.clear();
    for (Primary &prim : m_primaries) m_junctionEnds.push_back({ &prim.pipe, +1 });
    m_junctionEnds.push_back({ m_collector[0], -1 });
    m_junction.pressure = 0.0;

    m_outletMassFlow = 0.0;
    m_outletMassFlowDerivative = 0.0;
    m_massIn = m_massOut = 0.0;
    m_initialMass = totalMass();
    m_lastSubSteps = 0;
    m_stepOutflow = m_stepTime = 0.0;
    m_stepSubSteps = 0;
    m_built = true;
}

// Compressible orifice (reference: orifice), isentropic nozzle flow with the
// critical pressure ratio of the gas.
double ExhaustNetwork::orificeFlow(double CdA, double p0, double T0, double pd) const {
    if (!(CdA > 0) || !(p0 > pd) || !(T0 > 0)) return 0.0;
    const double g = m_gamma, g1 = g - 1;
    const double pr = pd / p0;
    const double rho0 = p0 / (m_R * T0);
    double psi;
    if (pr <= m_prc) {
        psi = m_psic;
    }
    else {
        // pr^(2/g) and pr^((g+1)/g) = pr^(2/g) * pr^((g-1)/g)
        const double a = std::pow(pr, 2 / g);
        const double b = a * std::pow(pr, g1 / g);
        psi = std::sqrt(std::max(0.0, 2 * g / g1 * (a - b)));
    }
    return CdA * std::sqrt(p0 * rho0) * psi;
}

// The port boundary at a primary's inlet face, solved before every pipe sub-step (the valve boundary of the
// method of characteristics: Benson, The Thermodynamics and Gas Dynamics of Internal-Combustion Engines, vol. 1,
// Oxford 1982; Winterbone & Pearson, Theory of Engine Manifold Design, 2000). The face pressure p_f is the
// one at which the flow the pipe takes along its incoming characteristic (u_f = J + 2 c_f / (gamma - 1), with
// J = u - 2 c / (gamma - 1) of the first cell and c_f isentropic from it) equals the compressible orifice flow
// between the port volume and the face:
//   port -> pipe (port pressure above the face's closed-end pressure): the gas at the face is port gas,
//     T_f = (h0_port - u_f^2 / 2) / cp, rho_f = p_f / (R T_f) (the contact sits at the face, as at the open end),
//     and the orifice runs from the port's (p, T) to p_f;
//   pipe -> port: the gas at the face is pipe gas, isentropic from the first cell to p_f, and the orifice runs
//     from its stagnation state to the port pressure.
// The port volume gives up or takes exactly the face's mass, enthalpy and species flux over the sub-step. For a
// small pressure difference the flow is then limited by the pipe's wave impedance (m ~ A dp / c) instead of the
// orifice's square root, whose slope is unbounded at dp = 0. The first version computed the orifice flow from
// the first cell's own state and added it to the cell as a source: between a 50 cc port and a 30-60 cc cell
// that explicit square-root exchange settles into a two-step oscillation of the pressure difference of about
// rho (C_d A dt c^2 (1/V_port + 1/V_cell))^2 / 2 (kPa at 20 kHz), so the flow reversed on nearly every
// sub-step and the results moved with the cell size and the simulation frequency.
void ExhaustNetwork::portBoundary(Primary &prim, double dts) {
    GasSystem *port = prim.spec.port;
    double F[Pipe1D::FluxCount];
    if (port == nullptr) {
        portFace(prim, nullptr, dts, 0.0, F, &prim.facePressure);
        prim.pipe.setEndFlux(-1, F);
        prim.lastMassFlow = 0.0;
        return;
    }
    // Explicit midpoint rule for the port volume (second order): the face flux from the port's state half a
    // sub-step on. The port relaxes toward the pipe on V / (c A), about 0.05 ms for a 50 cc port on a 44 mm pipe;
    // one explicit step per 6 us sub-step (20 kHz) makes it relax about 7% too fast, a first-order error that left
    // the torque drifting with the simulation frequency (the test I4 at 6000 rpm on 2 cm cells: 134.1 / 135.3 /
    // 140.6 N m at 15 / 20 / 40 kHz with one step, 137.8 / 138.5 / 139.5 with this rule).
    double F1[Pipe1D::FluxCount], pf1;
    portFace(prim, port, dts, prim.facePressure, F1, &pf1);
    GasSystem half = *port;
    applyPortFlux(half, F1, 0.5 * dts);
    portFace(prim, &half, dts, pf1, F, &prim.facePressure);
    prim.pipe.setEndFlux(-1, F);
    prim.lastMassFlow = F[0];
    applyPortFlux(*port, F, dts);
    m_massIn += F[0] * dts;
}

// Face flux of a primary's inlet for a given state of its port volume (see portBoundary); `port` null or empty, or
// no orifice: a closed end at the characteristic closed-end pressure. `guess` (a face pressure) starts the root
// search.
void ExhaustNetwork::portFace(
    const Primary &prim, const GasSystem *port, double dts, double guess, double *F, double *facePressure) const
{
    const Pipe1D &P = prim.pipe;
    const double g = m_gamma, g1 = g - 1, k2 = 2 / g1, ex = g1 / (2 * g), ig = 1 / g;
    const double A = P.endFaceArea(-1);
    const double CdA = std::min(prim.spec.orificeArea, A);   // never wider than the pipe it feeds

    const double rc = P.density(0), uc = P.velocity(0), pc = P.pressure(0);
    const double cc = std::sqrt(g * pc / rc);
    const double J = uc - k2 * cc;                 // left-running invariant reaching the face from the pipe
    const double cWall = -J / k2;                  // face sound speed (cell entropy) for u_f = 0
    const double pWall = cWall > 0 ? pc * std::pow(cWall / cc, 1 / ex) : 0.0;

    for (int k = 0; k < Pipe1D::FluxCount; ++k) F[k] = 0.0;
    const bool open = port != nullptr && CdA > 0 && port->n() > 0 && port->temperature() > 0;
    const double pr = open ? port->pressure() : 0.0;
    if (!open || pr == pWall) {
        F[1] = pWall * A;
        *facePressure = pWall;
        return;
    }
    const double Tr = port->temperature();
    const double M = m_params.molarMass;

    // h(pf) rises through zero between lo and hi. Secant steps from the guess (the last root), which is usually
    // within a fraction of a percent of the answer: three evaluations; Illinois (modified regula falsi) on the
    // bracket when a step leaves it or does not settle
    const double tolFlow = 1e-10 * CdA * pr / std::sqrt(m_R * Tr);
    auto solve = [&](auto &&h, double lo, double hi) {
        if (guess > lo && guess < hi) {
            double x0 = guess, h0 = h(x0);
            if (std::fabs(h0) < tolFlow) return x0;
            double x1 = x0 * (1 + 1e-6), h1 = h(x1);
            for (int it = 0; it < 4; ++it) {
                if (std::fabs(h1) < tolFlow) return x1;
                if (h1 == h0) break;
                const double x2 = x1 - h1 * (x1 - x0) / (h1 - h0);
                if (!(x2 > lo && x2 < hi)) break;
                x0 = x1; h0 = h1; x1 = x2; h1 = h(x1);
            }
        }
        double hlo = h(lo), hhi = h(hi);
        if (!(hlo < 0)) return lo;
        if (!(hhi > 0)) return hi;
        double x = 0.5 * (lo + hi);
        int side = 0;
        for (int it = 0; it < 60; ++it) {
            x = (lo * hhi - hi * hlo) / (hhi - hlo);
            if (!(x > lo && x < hi)) x = 0.5 * (lo + hi);
            const double hx = h(x);
            if (hx == 0 || hi - lo < 1e-11 * hi || std::fabs(hx) < tolFlow) return x;
            if (hx < 0) {
                lo = x; hlo = hx;
                if (side == -1) hhi *= 0.5;
                side = -1;
            }
            else {
                hi = x; hhi = hx;
                if (side == +1) hlo *= 0.5;
                side = +1;
            }
        }
        return x;
    };

    double pf, uf, f, h0;
    double x[Pipe1D::ScalarCount];
    if (pr > pWall) {
        // port -> pipe: port gas at the face; it carries its enthalpy and its share of the port's bulk kinetic
        // energy (the port loses that share with its momentum, so the pipe has to receive it)
        h0 = m_cp * Tr + (port->mass() > 0 ? port->bulkKineticEnergy() / port->mass() : 0.0);
        const double uStar = std::sqrt(2 * g1 * h0 / (g + 1));   // sonic at the face
        auto faceVelocity = [&](double p) { return std::min(uStar, J + k2 * cc * std::pow(p / pc, ex)); };
        auto pipeFlow = [&](double p) {
            const double u = faceVelocity(p);
            return p / (m_R * (h0 - 0.5 * u * u) / m_cp) * u * A;
        };
        pf = solve([&](double p) { return pipeFlow(p) - orificeFlow(CdA, pr, Tr, p); }, pWall, pr);
        uf = faceVelocity(pf);
        f = pipeFlow(pf);
        // at most 90% of the port's gas in one sub-step (as GasSystem::flow)
        f = std::min(f, MaxForwardMoleFraction * port->n() * M / dts);
        const GasSystem::Mix mix = port->mix();
        x[0] = mix.p_fuel;
        x[1] = mix.p_o2;
    }
    else {
        // pipe -> port: pipe gas at the face, isentropic from the first cell
        auto faceVelocity = [&](double p) {
            const double c = cc * std::pow(p / pc, ex);
            return std::max(-c, J + k2 * c);
        };
        auto state = [&](double p, double *T0, double *p0) {
            const double u = faceVelocity(p);
            const double T = p / (m_R * rc * std::pow(p / pc, ig));
            *T0 = T + 0.5 * u * u / m_cp;
            *p0 = p * std::pow(*T0 / T, g / g1);
            return -rc * std::pow(p / pc, ig) * u * A;      // kg/s out of the pipe
        };
        pf = solve([&](double p) {
            double T0, p0;
            const double out = state(p, &T0, &p0);
            return orificeFlow(CdA, p0, T0, pr) - out;
        }, pr, pWall);
        double T0, p0;
        f = -state(pf, &T0, &p0);
        uf = faceVelocity(pf);
        h0 = m_cp * T0;
        for (int k = 0; k < Pipe1D::ScalarCount; ++k) x[k] = P.scalar(0, k);
    }

    F[0] = f;
    F[1] = f * uf + pf * A;
    F[2] = f * h0;
    for (int k = 0; k < Pipe1D::ScalarCount; ++k) F[3 + k] = f * x[k];
    *facePressure = pf;
}

// The port volume takes the opposite of a face flux over dt: exactly -F[0] dt of mass, -F[2] dt of total energy
// (internal plus bulk kinetic) and the species with them.
void ExhaustNetwork::applyPortFlux(GasSystem &port, const double *F, double dt) const {
    const double M = m_params.molarMass;
    const double dm = F[0] * dt, dE = F[2] * dt;
    if (dm > 0 && port.n() > 0) {
        // gas leaves: it takes its share of the internal energy and of the bulk momentum (as GasSystem::flow; keeping
        // the momentum would speed up the gas left behind and make kinetic energy); whatever else the flux carries
        // (the flow work p v = R T per mole, the share of bulk kinetic energy) comes out of the internal energy
        const double n0 = port.n();
        const double dn = std::min(dm / M, n0);
        const double share = dn / n0;
        const double bulk = port.bulkKineticEnergy();
        const double perMol = port.kineticEnergyPerMol();
        // the species leave as the flux says (the port's own mix, unless the flux was taken half a sub-step on after
        // gas came in), so the port keeps the rest exactly
        GasSystem::Mix mix = port.mix();
        if (n0 - dn > 0) {
            mix.p_fuel = std::max(0.0, (n0 * mix.p_fuel - dn * F[3] / F[0]) / (n0 - dn));
            mix.p_o2 = std::max(0.0, (n0 * mix.p_o2 - dn * F[4] / F[0]) / (n0 - dn));
            mix.p_inert = std::max(0.0, 1.0 - mix.p_fuel - mix.p_o2);
        }
        port.scaleMomentum(1.0 - share);
        port.loseN(dn, perMol);
        port.changeEnergy(-(dE - dn * perMol - share * bulk));
        port.changeMix(mix);
    }
    else if (dm < 0) {
        // gas arrives with its stagnation enthalpy and no momentum along the port: the port's bulk motion slows and
        // the kinetic energy that loses becomes heat (as GasSystem::flow does for its sink)
        const double dn = -dm / M;
        GasSystem::Mix mix;
        mix.p_fuel = F[3] / F[0];
        mix.p_o2 = F[4] / F[0];
        mix.p_inert = std::max(0.0, 1.0 - mix.p_fuel - mix.p_o2);
        const double bulk0 = port.bulkKineticEnergy();
        port.gainN(dn, -dE / dn, mix);
        port.changeEnergy(bulk0 - port.bulkKineticEnergy());
    }
}

// Two pipes of different area meeting end to end (reference: junction).
// Mass and energy cross unchanged; the momentum fluxes differ by the force on
// the step.  Newton on the mass-flow mismatch, bisection as a fallback.
void ExhaustNetwork::junction2(Pipe1D &a, Pipe1D &b, Junction &jn) const {
    // the general junction with two ends: it also hands the gas crossing into the other pipe its own (upstream)
    // state, which the earlier two-pipe solver did not
    JunctionEnd ends[2] = { { &a, +1 }, { &b, -1 } };
    junctionN(ends, 2, jn);
}

// N pipe ends at one point (reference: junctionN): one pressure, mass flows
// that sum to zero, inflows mixed and sent down the receiving branches.
void ExhaustNetwork::junctionN(JunctionEnd *ends, int count, Junction &jn) const {
    const double g = m_gamma, g1 = g - 1;
    const double k2 = 2 / g1, ex = g1 / (2 * g), ig = 1 / g;

    struct EndState {
        double r, u, p, c, A, J;
        double x[Pipe1D::ScalarCount];
        double cs, us, rs, m;
    };
    EndState E[32];
    if (count > 32) count = 32;

    double pMin = 0, pMax = 0, pSum = 0;
    for (int k = 0; k < count; ++k) {
        const Pipe1D &P = *ends[k].pipe;
        const int side = ends[k].side, i = P.endCell(side);
        E[k].r = P.density(i); E[k].u = P.velocity(i); E[k].p = P.pressure(i);
        E[k].c = std::sqrt(g * E[k].p / E[k].r);
        E[k].A = P.endFaceArea(side);
        E[k].J = side > 0 ? E[k].u + k2 * E[k].c : E[k].u - k2 * E[k].c;
        for (int s = 0; s < Pipe1D::ScalarCount; ++s) E[k].x[s] = P.scalar(i, s);
        pMin = (k == 0) ? E[k].p : std::min(pMin, E[k].p);
        pMax = (k == 0) ? E[k].p : std::max(pMax, E[k].p);
        pSum += E[k].p;
    }

    // mass flow into the junction from end k at junction pressure ps, and its
    // derivative: s = ps / p_k, c = c_k s^ex, rho = rho_k s^ig = rho_k (s^ex)^k2,
    // d(rho u)/dps = rho (ig u -/+ k2 ex c) / ps
    auto inflow = [&](int k, double ps, double *derivative) {
        const int side = ends[k].side;
        const double se = std::pow(ps / E[k].p, ex);
        const double c = E[k].c * se;
        const double r = E[k].r * std::pow(se, k2);
        const double u = side > 0 ? E[k].J - k2 * c : E[k].J + k2 * c;
        E[k].cs = c; E[k].us = u; E[k].rs = r;
        if (derivative) {
            *derivative = side * E[k].A * r * (ig * u + (side > 0 ? -1.0 : 1.0) * k2 * ex * c) / ps;
        }
        return side * r * u * E[k].A;
    };
    auto mismatch = [&](double ps, double *derivative) {
        double t = 0, d = 0;
        for (int k = 0; k < count; ++k) {
            double dk;
            t += inflow(k, ps, derivative ? &dk : nullptr);
            if (derivative) d += dk;
        }
        if (derivative) *derivative = d;
        return t;
    };

    double ps = jn.pressure > 0 ? jn.pressure : pSum / count;
    bool ok = false;
    for (int it = 0; it < 10; ++it) {
        double d;
        const double gv = mismatch(ps, &d);
        if (!(d < 0)) break;
        const double nx = ps - gv / d;
        if (!(nx > 0)) break;
        if (std::fabs(nx - ps) < 1e-8 * ps) { ps = nx; ok = true; break; }
        ps = nx;
    }
    if (!ok) {
        double lo = pMin * 0.2, hi = pMax * 5;
        for (int it = 0; it < 60; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (mismatch(mid, nullptr) > 0) lo = mid; else hi = mid;
        }
        ps = 0.5 * (lo + hi);
    }

    // A receiving end (gas leaving the junction into that pipe) takes the mixed inflow's state at the face, as
    // the open end takes the outside air: the velocity from the pipe's own characteristic at ps, the temperature
    // from the mixed stagnation enthalpy (T = (h0mix - u^2/2) / cp) and the density ps / (R T). Giving it the
    // receiving pipe's own density made the mass flux disagree with the energy flux and put a spurious pressure
    // step at every temperature front (6-10 kPa for a 900 K pulse meeting 300 K gas). The Newton solution above,
    // with each end's own density, is the start; the mass balance with the mixed state is then solved by the
    // secant method (it falls monotonically with ps), with bisection as the fallback.
    auto balance = [&](double pj) {
        double mIn = 0, hIn = 0;
        for (int k = 0; k < count; ++k) {
            inflow(k, pj, nullptr);
            E[k].m = ends[k].side * E[k].rs * E[k].us * E[k].A;
            if (E[k].m > 0) {
                mIn += E[k].m;
                hIn += E[k].m * (E[k].cs * E[k].cs / g1 + 0.5 * E[k].us * E[k].us);
            }
        }
        double t = 0;
        for (int k = 0; k < count; ++k) {
            if (E[k].m < 0 && mIn > 0) {
                const double Tm = (hIn / mIn - 0.5 * E[k].us * E[k].us) / m_cp;
                if (Tm > 1.0) {
                    E[k].rs = pj / (m_R * Tm);
                    E[k].m = ends[k].side * E[k].rs * E[k].us * E[k].A;
                }
            }
            t += E[k].m;
        }
        return t;
    };
    double scale = 0;
    for (int k = 0; k < count; ++k) scale += E[k].A * E[k].r * E[k].c;
    const double tol = 1e-9 * scale;
    double pa = ps, fa = balance(pa);
    if (std::fabs(fa) > tol) {
        double pb = pa * (1 + 1e-4), fb = balance(pb);
        bool done = false;
        for (int it = 0; it < 25; ++it) {
            if (fb == fa) break;
            const double pn = pb - fb * (pb - pa) / (fb - fa);
            if (!(pn > 0)) break;
            pa = pb; fa = fb; pb = pn; fb = balance(pb);
            if (std::fabs(fb) <= tol || std::fabs(pb - pa) < 1e-10 * pb) { done = true; break; }
        }
        if (done) ps = pb;
        else {
            double lo = pMin * 0.2, hi = pMax * 5;
            for (int it = 0; it < 64; ++it) {
                const double mid = 0.5 * (lo + hi);
                if (balance(mid) > 0) lo = mid; else hi = mid;
            }
            ps = 0.5 * (lo + hi);
        }
    }
    jn.pressure = ps;
    balance(ps);
    // the solve leaves a residual of order 1e-9 of the flux scale: scale the receiving ends to the inflow so the
    // mass (and with it the energy and species) fluxes sum to zero to round-off
    {
        double mIn = 0, mOut = 0;
        for (int k = 0; k < count; ++k) {
            if (E[k].m > 0) mIn += E[k].m; else mOut -= E[k].m;
        }
        if (mIn > 0 && mOut > 0) {
            const double scaleOut = mIn / mOut;
            for (int k = 0; k < count; ++k) {
                if (E[k].m < 0) E[k].m *= scaleOut;
            }
        }
    }

    double mIn = 0, hIn = 0, xIn[Pipe1D::ScalarCount] = { 0.0, 0.0 };
    for (int k = 0; k < count; ++k) {
        if (E[k].m > 0) {
            const double h0 = E[k].cs * E[k].cs / g1 + 0.5 * E[k].us * E[k].us;
            mIn += E[k].m;
            hIn += E[k].m * h0;
            for (int s = 0; s < Pipe1D::ScalarCount; ++s) xIn[s] += E[k].m * E[k].x[s];
        }
    }
    const double hMix = mIn > 0 ? hIn / mIn : 0;
    double xMix[Pipe1D::ScalarCount];
    for (int s = 0; s < Pipe1D::ScalarCount; ++s) xMix[s] = mIn > 0 ? xIn[s] / mIn : 0;

    double F[Pipe1D::FluxCount];
    for (int k = 0; k < count; ++k) {
        const int side = ends[k].side;
        const double h0 = E[k].m > 0 ? E[k].cs * E[k].cs / g1 + 0.5 * E[k].us * E[k].us : hMix;
        const double f = side * E[k].m;              // flux in +x at that end
        F[0] = f;
        F[1] = f * E[k].us + ps * E[k].A;
        F[2] = f * h0;
        for (int s = 0; s < Pipe1D::ScalarCount; ++s) F[3 + s] = f * (E[k].m > 0 ? E[k].x[s] : xMix[s]);
        ends[k].pipe->setEndFlux(side, F);
    }
}

void ExhaustNetwork::applyBoundaries(double dts) {
    // the ports (the face flux of each primary's inlet, and the port volume's change over the sub-step)
    for (Primary &prim : m_primaries) portBoundary(prim, dts);
    junctionN(m_junctionEnds.data(), static_cast<int>(m_junctionEnds.size()), m_junction);
    for (size_t q = 1; q < m_collector.size(); ++q) {
        junction2(*m_collector[q - 1], *m_collector[q], m_collectorJunctions[q]);
    }
    m_collector.back()->applyOpenEnd(
        +1, m_params.atmosphericPressure, m_params.atmosphericTemperature, m_atmScalars);
}

void ExhaustNetwork::advance(double dt) {
    if (!m_built || !(dt > 0)) return;

    double dxMin = m_collector[0]->cellSize();
    for (const Primary &prim : m_primaries) dxMin = std::min(dxMin, prim.pipe.cellSize());
    for (const Pipe1D *p : m_collector) dxMin = std::min(dxMin, p->cellSize());

    Pipe1D *outlet = m_collector.back();
    double outflow = 0.0;
    double remaining = dt;
    int subSteps = 0;
    // The primitive state is kept current by Pipe1D::advance(); only the
    // initial state after build() or a setCellState() needs computePrimitives().
    while (remaining > 1e-15) {
        double smax = 0.0;
        for (Primary &prim : m_primaries) {
            smax = std::max(smax, prim.pipe.maxSignalSpeed());
        }
        for (Pipe1D *p : m_collector) {
            smax = std::max(smax, p->maxSignalSpeed());
        }
        const double dtCfl = m_params.cfl * dxMin / std::max(smax, 1.0);
        const double dts = std::min(remaining, dtCfl);

        applyBoundaries(dts);
        for (Primary &prim : m_primaries) prim.pipe.advance(dts);
        for (Pipe1D *p : m_collector) p->advance(dts);

        outflow += outlet->endMassFlow(+1) * dts;
        remaining -= dts;
        ++subSteps;
        if (subSteps > 10000) break;   // safety against a stalled CFL loop
    }

    m_stepOutflow += outflow;
    m_stepTime += dt;
    m_stepSubSteps += subSteps;
}

void ExhaustNetwork::finishStep() {
    if (!m_built || !(m_stepTime > 0)) return;
    const double massFlow = m_stepOutflow / m_stepTime;
    m_outletMassFlowDerivative = (massFlow - m_outletMassFlow) / m_stepTime;
    m_outletMassFlow = massFlow;
    m_collector.back()->setAveragedEndMassFlow(+1, massFlow);
    m_massOut += m_stepOutflow;
    m_lastSubSteps = m_stepSubSteps;
    m_stepOutflow = m_stepTime = 0.0;
    m_stepSubSteps = 0;
}

void ExhaustNetwork::step(double dt) {
    advance(dt);
    finishStep();
}

double ExhaustNetwork::outletPressure() const {
    const Pipe1D *p = m_collector.back();
    return p->pressure(p->cellCount() - 1);
}

double ExhaustNetwork::outletTemperature() const {
    const Pipe1D *p = m_collector.back();
    return p->temperature(p->cellCount() - 1);
}

GasSystem::Mix ExhaustNetwork::outletMix() const {
    const Pipe1D *p = m_collector.back();
    const int i = p->cellCount() - 1;
    GasSystem::Mix mix;
    mix.p_fuel = p->scalar(i, 0);
    mix.p_o2 = p->scalar(i, 1);
    mix.p_inert = std::max(0.0, 1.0 - mix.p_fuel - mix.p_o2);
    return mix;
}

double ExhaustNetwork::totalMass() const {
    double m = 0.0;
    for (const Primary &prim : m_primaries) m += prim.pipe.totalMass();
    for (const Pipe1D *p : m_collector) m += p->totalMass();
    return m;
}

double ExhaustNetwork::totalEnergy() const {
    double e = 0.0;
    for (const Primary &prim : m_primaries) e += prim.pipe.totalEnergy();
    for (const Pipe1D *p : m_collector) e += p->totalEnergy();
    return e;
}

int ExhaustNetwork::cellCount() const {
    int n = 0;
    for (const Primary &prim : m_primaries) n += prim.pipe.cellCount();
    for (const Pipe1D *p : m_collector) n += p->cellCount();
    return n;
}
