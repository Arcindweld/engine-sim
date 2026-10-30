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
    // Largest share of the end cell's mass that one port exchange may draw back
    // into the port volume; keeps the explicit coupling bounded when the pipe
    // end is at a much higher pressure than the port.
    constexpr double MaxReverseCellFraction = 0.05;
    // Same role as the 0.9 clamp in GasSystem::flow for the port volume.
    constexpr double MaxForwardMoleFraction = 0.9;
}

ExhaustNetwork::ExhaustNetwork() {
    m_built = false;
    m_gamma = 1.4;
    m_R = 287.0;
    m_cp = m_gamma / (m_gamma - 1) * m_R;
    m_atmScalars[0] = m_atmScalars[1] = 0.0;
    m_outletMassFlow = 0.0;
    m_outletMassFlowDerivative = 0.0;
    m_massIn = m_massOut = m_initialMass = 0.0;
    m_lastSubSteps = 0;
}

ExhaustNetwork::~ExhaustNetwork() {
    for (Pipe1D *p : m_collector) delete p;
    m_collector.clear();
}

int ExhaustNetwork::addPrimary(const PrimarySpec &spec) {
    Primary prim;
    prim.spec = spec;
    for (int k = 0; k < Pipe1D::FluxCount; ++k) prim.accumulated[k] = prim.averaged[k] = 0.0;
    prim.jetSpeed = 0.0;
    prim.lastMassFlow = 0.0;
    m_primaries.push_back(prim);
    return static_cast<int>(m_primaries.size()) - 1;
}

void ExhaustNetwork::build(const Parameters &params) {
    m_params = params;
    m_gamma = 1.0 + 2.0 / params.degreesOfFreedom;
    m_R = constants::R / params.molarMass;
    m_cp = m_gamma / (m_gamma - 1) * m_R;

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
        s.length = prim.spec.length;
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
    m_built = true;
}

// Compressible orifice (reference: orifice), isentropic nozzle flow with the
// critical pressure ratio of the gas.
double ExhaustNetwork::orificeFlow(double CdA, double p0, double T0, double pd) const {
    if (!(CdA > 0) || !(p0 > pd) || !(T0 > 0)) return 0.0;
    const double g = m_gamma, g1 = g - 1;
    const double prc = std::pow(2 / (g + 1), g / g1);
    const double pr = pd / p0;
    const double rho0 = p0 / (m_R * T0);
    double psi;
    if (pr <= prc) {
        psi = std::sqrt(g) * std::pow(2 / (g + 1), (g + 1) / (2 * g1));
    }
    else {
        psi = std::sqrt(std::max(0.0, 2 * g / g1 * (std::pow(pr, 2 / g) - std::pow(pr, (g + 1) / g))));
    }
    return CdA * std::sqrt(p0 * rho0) * psi;
}

double ExhaustNetwork::jetSpeed(double p0, double T0, double pd) const {
    const double g = m_gamma, g1 = g - 1;
    const double prc = std::pow(2 / (g + 1), g / g1);
    return std::sqrt(std::max(0.0, 2 * m_cp * T0 * (1 - std::pow(std::max(pd / p0, prc), g1 / g))));
}

// Flow between the port GasSystem and the first cell of the primary
// (reference: valveFlow + fedEnd).  Positive mass flow = into the pipe.
void ExhaustNetwork::exchangePort(int index, double dt) {
    Primary &prim = m_primaries[index];
    GasSystem *port = prim.spec.port;
    Pipe1D &P = prim.pipe;
    const int i = 0;
    const double g = m_gamma, g1 = g - 1;

    const double r = P.density(i), u = P.velocity(i), p = P.pressure(i);
    const double T = p / (m_R * r);
    const double T0 = T + 0.5 * u * u / m_cp;
    const double p0 = p * std::pow(T0 / T, g / g1);

    const double pr = port->pressure();
    const double Tr = port->temperature();
    const double M = m_params.molarMass;
    const double CdA = prim.spec.orificeArea;
    const double af = P.endFaceArea(-1);

    double m = 0.0, h = 0.0, uj = 0.0;
    double x[Pipe1D::ScalarCount] = { 0.0, 0.0 };
    GasSystem::Mix mix;

    if (p0 > pr) {
        // pipe -> port
        m = -orificeFlow(CdA, p0, T0, pr);
        const double cellMass = r * P.cellArea(i) * P.cellSize();
        const double limit = MaxReverseCellFraction * cellMass / dt;
        if (-m > limit) m = -limit;
        h = m_cp * T0;
        for (int k = 0; k < Pipe1D::ScalarCount; ++k) x[k] = P.scalar(i, k);
    }
    else if (port->n() > 0 && Tr > 0) {
        // port -> pipe
        m = orificeFlow(CdA, pr, Tr, p);
        const double limit = MaxForwardMoleFraction * port->n() * M / dt;
        if (m > limit) m = limit;
        h = m_cp * Tr;
        uj = jetSpeed(pr, Tr, p);
        mix = port->mix();
        x[0] = mix.p_fuel;
        x[1] = mix.p_o2;
    }

    if (m > 0) {
        const double dn = m * dt / M;
        // internal energy leaves with the gas (as in GasSystem::flow) plus the
        // flow work p v = R T per mole, so the pipe receives the full enthalpy
        port->loseN(dn, port->kineticEnergyPerMol());
        port->changeEnergy(-dn * constants::R * Tr);
    }
    else if (m < 0) {
        const double dn = -m * dt / M;
        mix.p_fuel = x[0];
        mix.p_o2 = x[1];
        mix.p_inert = std::max(0.0, 1.0 - x[0] - x[1]);
        // the port is a fixed volume: the incoming enthalpy becomes internal energy
        port->gainN(dn, h * M, mix);
    }

    const bool entering = m > 0;
    const double uf = entering ? uj : m / (r * af);
    prim.accumulated[0] += m * dt;
    prim.accumulated[1] += m * uf * dt;
    prim.accumulated[2] += m * h * dt;
    for (int k = 0; k < Pipe1D::ScalarCount; ++k) prim.accumulated[3 + k] += m * x[k] * dt;
    if (entering) prim.jetSpeed = std::max(prim.jetSpeed, uj);
    prim.lastMassFlow = m;
    m_massIn += m * dt;
}

void ExhaustNetwork::applyPortFlux(Primary &prim) {
    double F[Pipe1D::FluxCount];
    F[0] = prim.averaged[0];
    F[1] = prim.averaged[1] + prim.pipe.pressure(0) * prim.pipe.endFaceArea(-1);
    F[2] = prim.averaged[2];
    for (int k = 0; k < Pipe1D::ScalarCount; ++k) F[3 + k] = prim.averaged[3 + k];
    prim.pipe.setEndFlux(-1, F);
}

// Two pipes of different area meeting end to end (reference: junction).
// Mass and energy cross unchanged; the momentum fluxes differ by the force on
// the step.  Newton on the mass-flow mismatch, bisection as a fallback.
void ExhaustNetwork::junction2(Pipe1D &a, Pipe1D &b, Junction &jn) const {
    const int i = a.cellCount() - 1;
    const double A1 = a.endFaceArea(+1), A2 = b.endFaceArea(-1);
    const double rL = a.density(i), uL = a.velocity(i), pL = a.pressure(i);
    const double rR = b.density(0), uR = b.velocity(0), pR = b.pressure(0);
    const double g = m_gamma, g1 = g - 1;
    const double cL = std::sqrt(g * pL / rL), cR = std::sqrt(g * pR / rR);
    const double k2 = 2 / g1, ex = g1 / (2 * g), ig = 1 / g;
    const double JL = uL + k2 * cL, JR = uR - k2 * cR;

    auto mismatch = [&](double ps) {
        const double sL = ps / pL, sR = ps / pR;
        const double c1 = cL * std::pow(sL, ex), c2 = cR * std::pow(sR, ex);
        return A1 * rL * std::pow(sL, ig) * (JL - k2 * c1) - A2 * rR * std::pow(sR, ig) * (JR + k2 * c2);
    };

    double ps = jn.pressure > 0 ? jn.pressure : 0.5 * (pL + pR);
    bool ok = false;
    for (int it = 0; it < 8; ++it) {
        const double gv = mismatch(ps), hh = ps * 1e-6;
        const double d = (mismatch(ps + hh) - gv) / hh;
        if (!(d < 0)) break;
        const double nx = ps - gv / d;
        if (!(nx > 0)) break;
        if (std::fabs(nx - ps) < 1e-7 * ps) { ps = nx; ok = true; break; }
        ps = nx;
    }
    if (!ok || !(std::fabs(mismatch(ps)) <= 1e-6 * std::fabs(A1 * rL * cL))) {
        double lo = std::min(pL, pR) * 0.2, hi = std::max(pL, pR) * 5;
        for (int it = 0; it < 60; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (mismatch(mid) > 0) lo = mid; else hi = mid;
        }
        ps = 0.5 * (lo + hi);
    }
    jn.pressure = ps;

    const double c1 = cL * std::pow(ps / pL, ex), u1 = JL - k2 * c1, r1 = rL * std::pow(ps / pL, ig);
    const double c2 = cR * std::pow(ps / pR, ex), u2 = JR + k2 * c2;
    const double m = A1 * r1 * u1;                       // kg/s, + toward b
    const bool up = m >= 0;
    const double h0 = up ? c1 * c1 / g1 + 0.5 * u1 * u1 : c2 * c2 / g1 + 0.5 * u2 * u2;

    double Fa[Pipe1D::FluxCount], Fb[Pipe1D::FluxCount];
    Fa[0] = m; Fa[1] = m * u1 + ps * A1; Fa[2] = m * h0;
    Fb[0] = m; Fb[1] = m * u2 + ps * A2; Fb[2] = m * h0;
    for (int k = 0; k < Pipe1D::ScalarCount; ++k) {
        const double xx = up ? a.scalar(i, k) : b.scalar(0, k);
        Fa[3 + k] = Fb[3 + k] = m * xx;
    }
    a.setEndFlux(+1, Fa);
    b.setEndFlux(-1, Fb);
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

    // mass flow into the junction from end k at junction pressure ps
    auto inflow = [&](int k, double ps) {
        const int side = ends[k].side;
        const double s = ps / E[k].p;
        const double c = E[k].c * std::pow(s, ex);
        const double r = E[k].r * std::pow(s, ig);
        const double u = side > 0 ? E[k].J - k2 * c : E[k].J + k2 * c;
        E[k].cs = c; E[k].us = u; E[k].rs = r;
        return side * r * u * E[k].A;
    };
    auto mismatch = [&](double ps) {
        double t = 0;
        for (int k = 0; k < count; ++k) t += inflow(k, ps);
        return t;
    };

    double ps = jn.pressure > 0 ? jn.pressure : pSum / count;
    bool ok = false;
    for (int it = 0; it < 10; ++it) {
        const double gv = mismatch(ps), hh = ps * 1e-6;
        const double d = (mismatch(ps + hh) - gv) / hh;
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
            if (mismatch(mid) > 0) lo = mid; else hi = mid;
        }
        ps = 0.5 * (lo + hi);
    }
    jn.pressure = ps;
    mismatch(ps);

    double mIn = 0, hIn = 0, xIn[Pipe1D::ScalarCount] = { 0.0, 0.0 };
    for (int k = 0; k < count; ++k) {
        E[k].m = ends[k].side * E[k].rs * E[k].us * E[k].A;
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

void ExhaustNetwork::applyBoundaries() {
    for (Primary &prim : m_primaries) applyPortFlux(prim);
    junctionN(m_junctionEnds.data(), static_cast<int>(m_junctionEnds.size()), m_junction);
    for (size_t q = 1; q < m_collector.size(); ++q) {
        junction2(*m_collector[q - 1], *m_collector[q], m_collectorJunctions[q]);
    }
    m_collector.back()->applyOpenEnd(
        +1, m_params.atmosphericPressure, m_params.atmosphericTemperature, m_atmScalars);
}

void ExhaustNetwork::step(double dt) {
    if (!m_built || !(dt > 0)) return;

    for (Primary &prim : m_primaries) {
        for (int k = 0; k < Pipe1D::FluxCount; ++k) {
            prim.averaged[k] = prim.accumulated[k] / dt;
            prim.accumulated[k] = 0.0;
        }
    }

    double dxMin = m_collector[0]->cellSize();
    for (const Primary &prim : m_primaries) dxMin = std::min(dxMin, prim.pipe.cellSize());
    for (const Pipe1D *p : m_collector) dxMin = std::min(dxMin, p->cellSize());

    Pipe1D *outlet = m_collector.back();
    double outflow = 0.0;
    double remaining = dt;
    int subSteps = 0;
    while (remaining > 1e-15) {
        double smax = 0.0;
        for (Primary &prim : m_primaries) {
            prim.pipe.computePrimitives();
            smax = std::max(smax, prim.pipe.maxSignalSpeed());
            const double c0 = prim.pipe.soundSpeed(prim.pipe.pressure(0), prim.pipe.density(0));
            smax = std::max(smax, prim.jetSpeed + c0);
        }
        for (Pipe1D *p : m_collector) {
            p->computePrimitives();
            smax = std::max(smax, p->maxSignalSpeed());
        }
        const double dtCfl = m_params.cfl * dxMin / std::max(smax, 1.0);
        const double dts = std::min(remaining, dtCfl);

        applyBoundaries();
        for (Primary &prim : m_primaries) prim.pipe.advance(dts);
        for (Pipe1D *p : m_collector) p->advance(dts);

        outflow += outlet->endMassFlow(+1) * dts;
        remaining -= dts;
        ++subSteps;
        if (subSteps > 10000) break;   // safety against a stalled CFL loop
    }
    for (Primary &prim : m_primaries) {
        prim.pipe.computePrimitives();
        prim.jetSpeed = 0.0;
    }
    for (Pipe1D *p : m_collector) p->computePrimitives();

    const double massFlow = outflow / dt;
    m_outletMassFlowDerivative = (massFlow - m_outletMassFlow) / dt;
    m_outletMassFlow = massFlow;
    outlet->setAveragedEndMassFlow(+1, massFlow);
    m_massOut += outflow;
    m_lastSubSteps = subSteps;
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
