// Tests for the 1-D exhaust gas-dynamics model.
//
// Copyright (c) 2026 Arcindweld. MIT License (see LICENSE).

#include <gtest/gtest.h>

#include "../include/pipe_1d.h"
#include "../include/exhaust_network.h"
#include "../include/gas_system.h"
#include "../include/units.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

constexpr double Patm = 101325.0;
constexpr double Tatm = 298.15;
constexpr double AtmScalars[Pipe1D::ScalarCount] = { 0.0, 0.0 };

Pipe1D::Parameters straightPipe(double length, double diameter, double cellSize) {
    Pipe1D::Parameters p;
    Pipe1D::Section s;
    s.length = length;
    s.inletDiameter = s.outletDiameter = diameter;
    p.sections.push_back(s);
    p.cellSize = cellSize;
    p.maxCells = 100000;
    p.wallTemperature = Tatm;        // no heat transfer source at rest
    p.initialPressure = Patm;
    p.initialTemperature = Tatm;
    return p;
}

double cflStep(const Pipe1D &pipe) {
    return 0.8 * pipe.cellSize() / pipe.maxSignalSpeed();
}

// A transmissive (non-reflecting for these tests) end: the face flux is the Euler flux of the end cell itself.
void transmissiveEnd(Pipe1D &pipe, int side) {
    const Pipe1D::State st = pipe.state(pipe.endCell(side));
    double F[Pipe1D::FluxCount];
    pipe.hllc(st, st, F);
    for (int k = 0; k < Pipe1D::FluxCount; ++k) F[k] *= pipe.endFaceArea(side);
    pipe.setEndFlux(side, F);
}

// Exact solution of the Riemann problem for an ideal gas (Toro, Riemann Solvers and Numerical Methods for Fluid
// Dynamics, 3rd ed., 2009, ch. 4): the pressure functions f_K and their derivatives, Newton for the star pressure,
// and the wave pattern sampled at S = x / t.
struct Prim { double r, u, p; };

double pressureFunction(double p, const Prim &K, double g, double *derivative) {
    const double c = std::sqrt(g * K.p / K.r);
    if (p > K.p) {   // shock
        const double A = 2 / ((g + 1) * K.r), B = (g - 1) / (g + 1) * K.p;
        const double q = std::sqrt(A / (p + B));
        *derivative = q * (1 - 0.5 * (p - K.p) / (B + p));
        return (p - K.p) * q;
    }
    const double pr = p / K.p;   // rarefaction
    *derivative = 1 / (K.r * c) * std::pow(pr, -(g + 1) / (2 * g));
    return 2 * c / (g - 1) * (std::pow(pr, (g - 1) / (2 * g)) - 1);
}

void starState(const Prim &L, const Prim &R, double g, double *pStar, double *uStar) {
    double p = 0.5 * (L.p + R.p), dL, dR;
    for (int it = 0; it < 100; ++it) {
        const double f = pressureFunction(p, L, g, &dL) + pressureFunction(p, R, g, &dR) + R.u - L.u;
        const double next = std::max(1e-12, p - f / (dL + dR));
        const bool done = std::fabs(next - p) < 1e-14 * p;
        p = next;
        if (done) break;
    }
    *pStar = p;
    *uStar = 0.5 * (L.u + R.u) + 0.5 * (pressureFunction(p, R, g, &dR) - pressureFunction(p, L, g, &dL));
}

Prim sampleExact(const Prim &L, const Prim &R, double g, double S) {
    double p, u;
    starState(L, R, g, &p, &u);
    const double cL = std::sqrt(g * L.p / L.r), cR = std::sqrt(g * R.p / R.r);
    const double gm = (g - 1) / (g + 1);
    if (S <= u) {
        if (p > L.p) {
            const double SL = L.u - cL * std::sqrt((g + 1) / (2 * g) * p / L.p + (g - 1) / (2 * g));
            if (S <= SL) return L;
            return { L.r * (p / L.p + gm) / (gm * p / L.p + 1), u, p };
        }
        const double head = L.u - cL, tail = u - cL * std::pow(p / L.p, (g - 1) / (2 * g));
        if (S <= head) return L;
        if (S >= tail) return { L.r * std::pow(p / L.p, 1 / g), u, p };
        const double k = 2 / (g + 1) + gm / cL * (L.u - S);
        return { L.r * std::pow(k, 2 / (g - 1)), 2 / (g + 1) * (cL + (g - 1) / 2 * L.u + S),
                 L.p * std::pow(k, 2 * g / (g - 1)) };
    }
    if (p > R.p) {
        const double SR = R.u + cR * std::sqrt((g + 1) / (2 * g) * p / R.p + (g - 1) / (2 * g));
        if (S >= SR) return R;
        return { R.r * (p / R.p + gm) / (gm * p / R.p + 1), u, p };
    }
    const double head = R.u + cR, tail = u + cR * std::pow(p / R.p, (g - 1) / (2 * g));
    if (S >= head) return R;
    if (S <= tail) return { R.r * std::pow(p / R.p, 1 / g), u, p };
    const double k = 2 / (g + 1) - gm / cR * (R.u - S);
    return { R.r * std::pow(k, 2 / (g - 1)), 2 / (g + 1) * (-cR + (g - 1) / 2 * R.u + S),
             R.p * std::pow(k, 2 * g / (g - 1)) };
}

struct ShockTube {
    double l1Rho, l1P, l1U;        // L1 errors relative to the L1 norm of the exact rho, p and max |u|
    double plateauP, pStar;        // mean computed pressure inside the star region, and the exact p*
    double massDrift;
};

// A 1 m tube (no friction, no heat transfer) split at 0.5 m, run to tEnd on n cells; the waves stay clear of the
// ends, which are transmissive.
ShockTube runShockTube(const Prim &L, const Prim &R, double tEnd, int n) {
    const double g = 1.4, x0 = 0.5;
    Pipe1D::Parameters pp = straightPipe(1.0, 0.05, 1.0 / n);
    pp.frictionFactor = 0.0;
    Pipe1D pipe;
    pipe.initialize(pp);
    for (int i = 0; i < pipe.cellCount(); ++i) {
        const Prim &w = (i + 0.5) * pipe.cellSize() < x0 ? L : R;
        pipe.setCellState(i, w.r, w.u, w.p, AtmScalars);
    }
    pipe.computePrimitives();
    const double m0 = pipe.totalMass();
    double t = 0;
    while (t < tEnd - 1e-15) {
        transmissiveEnd(pipe, -1);
        transmissiveEnd(pipe, +1);
        const double dt = std::min(cflStep(pipe), tEnd - t);
        pipe.advance(dt);
        t += dt;
    }
    ShockTube res{};
    double pStar, uStar, rhoRef = 0, pRef = 0, uRef = 1e-30;
    starState(L, R, g, &pStar, &uStar);
    std::vector<int> plateau;
    for (int i = 0; i < pipe.cellCount(); ++i) {
        const Prim e = sampleExact(L, R, g, ((i + 0.5) * pipe.cellSize() - x0) / t);
        res.l1Rho += std::fabs(pipe.density(i) - e.r);
        res.l1P += std::fabs(pipe.pressure(i) - e.p);
        res.l1U += std::fabs(pipe.velocity(i) - e.u);
        rhoRef += e.r; pRef += e.p; uRef = std::max(uRef, std::fabs(e.u));
        if (std::fabs(e.p - pStar) < 1e-9 * pStar) plateau.push_back(i);
    }
    res.l1Rho /= rhoRef;
    res.l1P /= pRef;
    res.l1U /= uRef * pipe.cellCount();
    // the star region without the cells next to its edges (the captured waves are a few cells wide)
    double sum = 0;
    int count = 0;
    for (size_t k = 4; k + 4 < plateau.size(); ++k) { sum += pipe.pressure(plateau[k]); ++count; }
    res.plateauP = count > 0 ? sum / count : 0;
    res.pStar = pStar;
    res.massDrift = (pipe.totalMass() - m0) / m0;
    return res;
}

} // namespace

// (a) A pipe at rest stays at rest: closed at one end, open at the other,
// uniform atmospheric state.  The open end must not inject or drain anything.
TEST(Pipe1DTests, PipeAtRestStaysAtRest) {
    Pipe1D pipe;
    pipe.initialize(straightPipe(1.0, 0.05, 0.02));
    const double m0 = pipe.totalMass();
    const double e0 = pipe.totalEnergy();
    for (int step = 0; step < 2000; ++step) {
        pipe.computePrimitives();
        pipe.applyWallEnd(-1);
        pipe.applyOpenEnd(+1, Patm, Tatm, AtmScalars);
        pipe.advance(cflStep(pipe));
    }
    pipe.computePrimitives();
    for (int i = 0; i < pipe.cellCount(); ++i) {
        EXPECT_NEAR(pipe.pressure(i), Patm, 1e-6 * Patm);
        EXPECT_NEAR(pipe.velocity(i), 0.0, 1e-9);
        EXPECT_NEAR(pipe.temperature(i), Tatm, 1e-6 * Tatm);
    }
    EXPECT_NEAR(pipe.totalMass(), m0, 1e-9 * m0);
    EXPECT_NEAR(pipe.totalEnergy(), e0, 1e-9 * e0);
}

// (b) A small pressure step in a closed-open pipe: the pulse travels to the
// open end, reflects inverted and comes back to the closed end (where it
// reflects with the same sign).  Check the sign and the round-trip time 2L/c
// within 2%.
TEST(Pipe1DTests, ClosedOpenPipePulseReflection) {
    const double L = 1.0, dx = 0.01;
    Pipe1D pipe;
    pipe.initialize(straightPipe(L, 0.05, dx));
    const double c = std::sqrt(1.4 * 287.0 * Tatm);  // expected for the default gas model
    EXPECT_NEAR(pipe.soundSpeed(Patm, Patm / (287.0 * Tatm)), c, 1e-9 * c);

    // Smooth compression pulse of +2 kPa centred on the closed end, at rest.
    // Half of it runs toward the open end; the other half is its mirror image
    // in the wall, so the two coincide: one pulse travels L to the open end,
    // reflects inverted and arrives back at the wall after 2L/c.
    const double dp = 2000.0, width = 0.1, x0 = 0.0;
    const double rho0 = Patm / (287.0 * Tatm);
    for (int i = 0; i < pipe.cellCount(); ++i) {
        const double x = (i + 0.5) * dx;
        const double g = std::exp(-0.5 * (x - x0) * (x - x0) / (width * width / 4));
        const double p = Patm + dp * g;
        const double rho = rho0 * std::pow(p / Patm, 1 / 1.4);
        pipe.setCellState(i, rho, 0.0, p, AtmScalars);
    }

    double t = 0;
    double tMin = 0, pMin = 1e30, pMax = -1e30;
    const int probe = 0;   // closed end
    while (t < 2.6 * L / c) {
        pipe.computePrimitives();
        pipe.applyWallEnd(-1);
        pipe.applyOpenEnd(+1, Patm, Tatm, AtmScalars);
        const double dt = cflStep(pipe);
        pipe.advance(dt);
        t += dt;
        pipe.computePrimitives();
        const double p = pipe.pressure(probe);
        // ignore the initial pulse itself (first 0.4 L/c)
        if (t > 0.8 * L / c) {
            if (p < pMin) { pMin = p; tMin = t; }
            if (p > pMax) pMax = p;
        }
    }
    // Reflection from the open end inverts the pulse: a rarefaction arrives at the wall.
    EXPECT_LT(pMin - Patm, -0.3 * dp);          // open-end reflection coefficient ~ -1 (doubled at the wall)
    EXPECT_LT(pMax - Patm, 0.3 * dp);           // no comparable compression came back
    const double expected = 2 * L / c;
    std::cout << "[  values  ] reflected minimum " << (pMin - Patm) << " Pa at " << tMin * 1e3
              << " ms, expected 2L/c = " << expected * 1e3 << " ms (" << (tMin / expected - 1) * 100 << " %)\n";
    EXPECT_NEAR(tMin, expected, 0.02 * expected);
}

// (c) Quarter-wave resonance of a closed-open pipe: excite with a broadband
// pulse, record the pressure at the closed end and find the dominant
// frequency; it must be within 3% of c / (4 L), with the usual end correction
// of an unflanged open pipe (0.61 r) taken into account.
TEST(Pipe1DTests, QuarterWaveResonance) {
    const double L = 0.8, d = 0.04, dx = 0.01;
    Pipe1D pipe;
    pipe.initialize(straightPipe(L, d, dx));
    const double c = std::sqrt(1.4 * 287.0 * Tatm);
    const double rho0 = Patm / (287.0 * Tatm);

    for (int i = 0; i < pipe.cellCount(); ++i) {
        const double x = (i + 0.5) * dx;
        const double g = std::exp(-0.5 * (x - 0.1) * (x - 0.1) / (0.03 * 0.03));
        const double p = Patm + 1000.0 * g;
        pipe.setCellState(i, rho0 * std::pow(p / Patm, 1 / 1.4), 0.0, p, AtmScalars);
    }

    // fixed time step (half the CFL step of the gas at rest; the pulse is
    // small so the signal speed hardly changes) and one sample per step
    pipe.computePrimitives();
    const double dt = 0.5 * cflStep(pipe);
    const double fs = 1.0 / dt, duration = 0.5;
    std::vector<double> samples;
    double t = 0;
    while (t < duration) {
        pipe.applyWallEnd(-1);
        pipe.applyOpenEnd(+1, Patm, Tatm, AtmScalars);
        pipe.advance(dt);
        t += dt;
        samples.push_back(pipe.pressure(0) - Patm);
    }

    // DFT magnitude search between 20 Hz and 500 Hz (0.5 Hz steps), Hann window
    double bestF = 0, bestMag = -1;
    const int n = static_cast<int>(samples.size());
    for (double f = 20; f <= 500; f += 0.5) {
        double re = 0, im = 0;
        for (int k = 0; k < n; ++k) {
            const double w = 0.5 - 0.5 * std::cos(2 * 3.14159265358979323846 * k / (n - 1.0));
            const double ph = 2 * 3.14159265358979323846 * f * k / fs;
            re += w * samples[k] * std::cos(ph);
            im -= w * samples[k] * std::sin(ph);
        }
        const double mag = re * re + im * im;
        if (mag > bestMag) { bestMag = mag; bestF = f; }
    }
    // The 1-D open-end boundary has no end correction (p = p_atm at the exit
    // face), so the ideal c / 4L is the reference.
    const double expected = c / (4 * L);
    std::cout << "[  values  ] resonance " << bestF << " Hz, c/4L = " << expected << " Hz ("
              << (bestF / expected - 1) * 100 << " %)\n";
    EXPECT_NEAR(bestF, expected, 0.03 * expected);
}

// (e) Junction of three pipes: a pulse arrives from one branch.  The
// junction keeps one pressure for all ends and the network conserves mass to
// round-off; the transmitted amplitude follows linear acoustics, 2 A1 / (A1 +
// A2 + A3) (Blair, Design and Simulation of Four-Stroke Engines, sec. 2.10;
// the same check the JS reference passes within 0.3%).
TEST(Pipe1DTests, JunctionPressureAndMassBalance) {
    ExhaustNetwork net;
    ExhaustNetwork::Parameters np;
    Pipe1D::Section s;
    s.length = 0.6;
    s.inletDiameter = s.outletDiameter = 0.06;
    np.collectorSections.push_back(s);
    np.cellSize = 0.005;
    np.maxCellsPerPipe = 100000;
    np.wallTemperature = Tatm;
    np.initialTemperature = Tatm;

    GasSystem portA, portB;
    portA.initialize(Patm, 1e-4, Tatm);
    portB.initialize(Patm, 1e-4, Tatm);
    ExhaustNetwork::PrimarySpec spec;
    spec.length = 0.6;
    spec.diameter = 0.04;
    spec.orificeArea = 0.0;   // ports closed for this test
    spec.port = &portA;
    net.addPrimary(spec);
    spec.port = &portB;
    net.addPrimary(spec);
    net.build(np);

    const double rho0 = Patm / (287.0 * Tatm);
    Pipe1D &a = net.primaryPipe(0);
    const Pipe1D &b = net.primaryPipe(1);
    const Pipe1D &coll = net.collectorPipe(0);
    for (int i = 0; i < a.cellCount(); ++i) {
        const double x = (i + 0.5) * a.cellSize();
        const double g = std::exp(-0.5 * (x - 0.15) * (x - 0.15) / (0.04 * 0.04));
        const double p = Patm + 2000.0 * g;
        a.setCellState(i, rho0 * std::pow(p / Patm, 1 / 1.4), 0.0, p, AtmScalars);
    }
    a.computePrimitives();
    const double m0 = net.totalMass();

    // probes: incident half-pulse in A at x = 0.4 m (before the junction
    // reflection comes back), transmitted pulses in B and the collector 0.2 m
    // past the junction
    const int probeA = static_cast<int>(0.4 / a.cellSize());
    const int probeB = b.cellCount() - 1 - static_cast<int>(0.2 / b.cellSize());
    const int probeC = static_cast<int>(0.2 / coll.cellSize());
    double peakA = 0, peakB = 0, peakC = 0;
    double t = 0;
    const double dt = 2e-5;
    while (t < 0.0025) {
        net.step(dt);
        t += dt;
        if (t < 0.0012) peakA = std::max(peakA, a.pressure(probeA) - Patm);
        peakB = std::max(peakB, b.pressure(probeB) - Patm);
        peakC = std::max(peakC, coll.pressure(probeC) - Patm);
        EXPECT_NEAR(net.totalMass() + net.massOut(), m0, 1e-9 * m0);
    }
    EXPECT_GT(peakA, 500.0);
    const double A1 = 0.04 * 0.04, A3 = 0.06 * 0.06;   // proportional to the areas
    const double transmission = 2 * A1 / (A1 + A1 + A3);
    std::cout << "[  values  ] junction transmission B " << peakB / peakA << ", collector " << peakC / peakA
              << ", linear theory " << transmission << "\n";
    EXPECT_NEAR(peakB / peakA, transmission, 0.06);
    EXPECT_NEAR(peakC / peakA, transmission, 0.06);
    // the collector and the second primary see the same junction pressure
    EXPECT_NEAR(peakB, peakC, 0.03 * peakA);
}

// Orifice flow of air through 1 cm^2 from 2 bar, 300 K (isentropic nozzle, e.g. Anderson, Modern Compressible
// Flow, ch. 5): choked below the critical pressure ratio, m = A p0 sqrt(g / (R T0)) (2/(g+1))^((g+1)/(2(g-1))),
// and subsonic above it, m = A p0 sqrt(2 g / ((g-1) R T0) (pr^(2/g) - pr^((g+1)/g))). The two meet at the
// critical ratio, the flow falls to zero as the pressures equalise, and nothing flows against the pressure drop.
// The gas is the network's after build(): gamma = 1 + 2 / 5, R = R_u / M of air.
TEST(Pipe1DTests, OrificeChokedFlow) {
    ExhaustNetwork net;
    GasSystem port;
    port.initialize(Patm, 1e-4, Tatm);
    ExhaustNetwork::PrimarySpec spec;
    spec.port = &port;
    spec.length = 0.3;
    spec.diameter = 0.04;
    net.addPrimary(spec);
    net.build(ExhaustNetwork::Parameters());

    const double g = 1.4, R = 8.31446261815324 / 0.02897;
    const double p0 = 2e5, T0 = 300.0, A = 1e-4;
    const double choked = A * p0 * std::sqrt(g / (R * T0)) * std::pow(2 / (g + 1), (g + 1) / (2 * (g - 1)));
    auto subsonic = [&](double pr) {
        return A * p0 * std::sqrt(2 * g / ((g - 1) * R * T0) * (std::pow(pr, 2 / g) - std::pow(pr, (g + 1) / g)));
    };
    const double prc = std::pow(2 / (g + 1), g / (g - 1));   // 0.528
    EXPECT_NEAR(net.orificeFlow(A, p0, T0, Patm), choked, 1e-9 * choked);
    EXPECT_NEAR(net.orificeFlow(A, p0, T0, 0.2 * p0), choked, 1e-9 * choked);
    for (double pr : { 0.6, 0.8, 0.9, 0.99 }) {
        EXPECT_NEAR(net.orificeFlow(A, p0, T0, pr * p0), subsonic(pr), 1e-9 * choked) << "pr " << pr;
    }
    EXPECT_NEAR(subsonic(prc), choked, 1e-9 * choked);   // the formulas agree at the critical ratio
    EXPECT_NEAR(net.orificeFlow(A, p0, T0, prc * p0 * (1 + 1e-9)), choked, 1e-6 * choked);
    double last = choked;
    for (double pr = 0.53; pr < 1.0; pr += 0.01) {
        const double m = net.orificeFlow(A, p0, T0, pr * p0);
        EXPECT_LT(m, last) << "pr " << pr;
        last = m;
    }
    EXPECT_EQ(net.orificeFlow(A, p0, T0, p0), 0.0);
    EXPECT_EQ(net.orificeFlow(A, Patm, T0, p0), 0.0);
}

// (h) Sod's shock tube (Sod, J. Comput. Phys. 27, 1978; here 1 / 0.125 kg/m3 and 1e5 / 1e4 Pa) and an exhaust
// blowdown (5 bar 1200 K into 1 bar 300 K) against the exact Riemann solution: the L1 errors are small on 400
// cells and fall as the grid is refined, the pressure between the rarefaction and the shock is the exact p*,
// and no mass is lost.
TEST(Pipe1DTests, ShockTubesMatchTheExactSolution) {
    struct Case { const char *name; Prim L, R; double tEnd, maxL1Rho, maxL1P; };
    const double Rg = 287.0;
    const Case cases[] = {
        { "Sod", { 1.0, 0, 1e5 }, { 0.125, 0, 1e4 }, 0.5e-3, 0.006, 0.005 },
        { "blowdown", { 5e5 / (Rg * 1200), 0, 5e5 }, { 1e5 / (Rg * 300), 0, 1e5 }, 0.4e-3, 0.015, 0.005 },
    };
    for (const Case &c : cases) {
        const ShockTube coarse = runShockTube(c.L, c.R, c.tEnd, 100);
        const ShockTube fine = runShockTube(c.L, c.R, c.tEnd, 400);
        std::cout << "[  values  ] " << c.name << ": L1 rho " << 100 * coarse.l1Rho << " % (100 cells) -> "
                  << 100 * fine.l1Rho << " % (400), p " << 100 * fine.l1P << " %, u " << 100 * fine.l1U
                  << " %; star p " << fine.plateauP << " Pa, exact " << fine.pStar << " Pa\n";
        EXPECT_LT(fine.l1Rho, c.maxL1Rho) << c.name;
        EXPECT_LT(fine.l1P, c.maxL1P) << c.name;
        EXPECT_LT(fine.l1Rho, 0.5 * coarse.l1Rho) << c.name;
        EXPECT_LT(fine.l1P, 0.5 * coarse.l1P) << c.name;
        EXPECT_NEAR(fine.plateauP, fine.pStar, 0.002 * fine.pStar) << c.name;
        EXPECT_NEAR(fine.massDrift, 0.0, 1e-12) << c.name;
    }
}

// (i) A temperature front (900 K gas meeting 300 K gas at the same pressure and velocity) carried through a
// junction. The exact answer is that only the temperature moves: pressure and velocity stay uniform. The
// junction hands each receiving pipe the arriving gas's own state; giving the receiving end its own density (the
// first version of the solver) put 6.5 kPa (flow one way) and 10 kPa (the other) pressure steps here, about half
// of rho c u. Two pipes through junction2, the same through a two-ended junctionN, and two primaries meeting a
// collector of their combined area (the mixing path), each with the flow both ways.
TEST(Pipe1DTests, TemperatureFrontCrossesJunctionsUndisturbed) {
    const double dx = 0.01, d = 0.04, tEnd = 1.0e-3, Rg = 287.0, Thot = 900.0, Tcold = 300.0;
    auto fill = [&](Pipe1D &p, double T, double u) {
        for (int i = 0; i < p.cellCount(); ++i) p.setCellState(i, Patm / (Rg * T), u, Patm, AtmScalars);
        p.computePrimitives();
    };
    auto makePipe = [&](Pipe1D &p, double diameter) {
        Pipe1D::Parameters pp = straightPipe(1.0, diameter, dx);
        pp.frictionFactor = 0.0;
        p.initialize(pp);
    };
    ExhaustNetwork net;   // only its junction solvers are used here (gamma 1.4, R 287 before build())
    for (int layout = 0; layout < 3; ++layout) {
        for (double u : { 50.0, -50.0 }) {
            Pipe1D a, b, c;
            makePipe(a, d);
            makePipe(b, d);
            makePipe(c, layout == 2 ? d * std::sqrt(2.0) : d);
            fill(a, Thot, u);
            fill(b, layout == 2 ? Thot : Tcold, u);
            fill(c, Tcold, u);
            ExhaustNetwork::Junction jn;
            ExhaustNetwork::JunctionEnd two[2] = { { &a, +1 }, { &c, -1 } };
            ExhaustNetwork::JunctionEnd three[3] = { { &a, +1 }, { &b, +1 }, { &c, -1 } };
            std::vector<Pipe1D *> pipes = { &a, &c };
            if (layout == 2) pipes.push_back(&b);
            double t = 0, dpMax = 0, duMax = 0;
            while (t < tEnd) {
                for (Pipe1D *q : pipes) transmissiveEnd(*q, q == &c ? +1 : -1);
                if (layout == 0) net.junction2(a, c, jn);
                else if (layout == 1) net.junctionN(two, 2, jn);
                else net.junctionN(three, 3, jn);
                double s = 0;
                for (Pipe1D *q : pipes) s = std::max(s, q->maxSignalSpeed());
                const double dt = std::min(0.8 * dx / s, tEnd - t);
                for (Pipe1D *q : pipes) q->advance(dt);
                t += dt;
                for (Pipe1D *q : pipes) {
                    for (int i = 0; i < q->cellCount(); ++i) {
                        dpMax = std::max(dpMax, std::fabs(q->pressure(i) - Patm));
                        duMax = std::max(duMax, std::fabs(q->velocity(i) - u));
                    }
                }
            }
            const char *names[3] = { "junction2", "junctionN, 2 ends", "junctionN, 2 primaries + collector" };
            std::cout << "[  values  ] " << names[layout] << ", u " << u << " m/s: max |p - p0| " << dpMax
                      << " Pa, max |u - u0| " << duMax << " m/s\n";
            EXPECT_LT(dpMax, 1.0) << names[layout] << ", u " << u;
            EXPECT_LT(duMax, 0.01) << names[layout] << ", u " << u;
            // the front really crossed: the receiving pipe's first cells changed temperature
            const Pipe1D &recv = u > 0 ? c : a;
            const int cell = u > 0 ? 0 : a.cellCount() - 1;
            const double Tarriving = u > 0 ? Thot : Tcold;
            EXPECT_NEAR(recv.temperature(cell), Tarriving, 0.05 * Tarriving) << names[layout] << ", u " << u;
        }
    }
}
