// Tests for the 1-D exhaust gas-dynamics model.
//
// Copyright (c) 2026 Arcindweld. MIT License (see LICENSE).

#include <gtest/gtest.h>

#include "../include/pipe_1d.h"
#include "../include/exhaust_network.h"
#include "../include/gas_system.h"
#include "../include/units.h"

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

// Orifice flow sanity: choked mass flow of air through 1 cm^2 at 2 bar, 300 K
// m = Cd A p0 sqrt(gamma / (R T0)) (2/(gamma+1))^((gamma+1)/(2(gamma-1)))
TEST(Pipe1DTests, OrificeChokedFlow) {
    ExhaustNetwork net;
    const double g = 1.4, R = 8.31446261815324 / 0.02897;
    const double p0 = 2e5, T0 = 300.0, A = 1e-4;
    const double expected = A * p0 * std::sqrt(g / (R * T0)) * std::pow(2 / (g + 1), (g + 1) / (2 * (g - 1)));
    EXPECT_NEAR(net.orificeFlow(A, p0, T0, Patm), expected, 1e-4 * expected);
    EXPECT_NEAR(net.orificeFlow(A, p0, T0, 0.9 * p0), net.orificeFlow(A, p0, T0, 0.9 * p0), 0.0);
    EXPECT_EQ(net.orificeFlow(A, Patm, T0, p0), 0.0);
}
