// Network-level tests for the 1-D exhaust model: conservation through the
// port coupling and the whole network.
//
// Copyright (c) 2026 Arcindweld. MIT License (see LICENSE).

#include <gtest/gtest.h>

#include "../include/exhaust_network.h"
#include "../include/gas_system.h"
#include "../include/units.h"

#include <cmath>
#include <iostream>

namespace {

constexpr double Patm = 101325.0;
constexpr double Tatm = 298.15;

// Four ports pulsed like a 4-cylinder at 6000 rpm: each port is refilled with
// hot gas during its "blowdown" window and left alone otherwise.
struct PulsedPorts {
    GasSystem port[4];
    double volume = 60e-6;

    PulsedPorts() {
        for (int i = 0; i < 4; ++i) {
            port[i].initialize(Patm, volume, Tatm);
            port[i].setGeometry(0.1, 0.03, 1.0, 0.0);
        }
    }

    // returns the mass (kg) added to the ports by this call
    double drive(double t, double dt) {
        const double period = 2.0 / (6000.0 / 60.0);   // one engine cycle, s
        double added = 0.0;
        for (int i = 0; i < 4; ++i) {
            const double phase = std::fmod(t + i * period / 4, period) / period;
            if (phase < 0.12) {
                // feed hot gas at a fixed molar rate (like a valve opening)
                const double nDot = 0.5;                     // mol/s
                GasSystem::Mix mix;
                mix.p_fuel = 0.0; mix.p_o2 = 0.01; mix.p_inert = 0.99;
                const double E = GasSystem::kineticEnergyPerMol(1200.0, 5) + 8.31446261815324 * 1200.0;
                port[i].gainN(nDot * dt, E, mix);
                added += nDot * dt * units::AirMolecularMass;
            }
        }
        return added;
    }
};

ExhaustNetwork::Parameters fourIntoOne() {
    ExhaustNetwork::Parameters np;
    Pipe1D::Section s;
    s.length = 0.6;
    s.inletDiameter = s.outletDiameter = 0.06;
    np.collectorSections.push_back(s);
    np.cellSize = 0.04;
    return np;
}

} // namespace

// (d) Mass conservation over 1 s of pulsed running: mass added to the ports
// = mass left in the ports + mass in the pipes + mass out of the open end,
// to far better than 0.5%.
TEST(ExhaustNetworkTests, MassConservationOneSecond) {
    PulsedPorts ports;
    ExhaustNetwork net;
    ExhaustNetwork::PrimarySpec spec;
    spec.length = 0.7;
    spec.diameter = 0.04;
    spec.orificeArea = 8e-4;
    for (int i = 0; i < 4; ++i) {
        spec.port = &ports.port[i];
        net.addPrimary(spec);
    }
    net.build(fourIntoOne());

    double portMass0 = 0.0;
    for (int i = 0; i < 4; ++i) portMass0 += ports.port[i].mass();
    const double pipeMass0 = net.totalMass();

    const double dt = 1e-4, dtFluid = dt / 8;
    double added = 0.0, t = 0.0;
    int steps = 0;
    while (t < 1.0) {
        for (int s = 0; s < 8; ++s) {
            added += ports.drive(t + s * dtFluid, dtFluid);
            for (int i = 0; i < 4; ++i) net.exchangePort(i, dtFluid);
        }
        net.step(dt);
        t += dt;
        ++steps;
        ASSERT_FALSE(std::isnan(net.totalMass())) << "NaN at t = " << t;
    }

    double portMass = 0.0;
    for (int i = 0; i < 4; ++i) portMass += ports.port[i].mass();
    const double pipeMass = net.totalMass();
    const double balance = (portMass - portMass0) + (pipeMass - pipeMass0) + net.massOut() - added;
    std::cout << "[  values  ] added " << added << " kg, out " << net.massOut() << " kg, pipes " << pipeMass
              << " kg, ports " << portMass << " kg, imbalance " << balance
              << " kg (" << 100 * balance / added << " % of the throughput), sub-steps/step "
              << net.lastSubSteps() << "\n";
    EXPECT_GT(added, 0.005);
    EXPECT_GT(net.massOut(), 0.5 * added);          // most of it left through the tailpipe
    EXPECT_NEAR(balance / added, 0.0, 0.005);
    EXPECT_NEAR(net.massIn(), added - (portMass - portMass0), 1e-6 * added);
}

// Energy conservation of the port exchange: a hot port emptying into a cold
// pipe.  Total (port internal energy + pipe energy) is conserved because the
// port gives up enthalpy and the pipe receives enthalpy + jet kinetic energy
// as total energy.
TEST(ExhaustNetworkTests, PortExchangeConservesEnergy) {
    GasSystem port;
    port.initialize(3.0 * Patm, 100e-6, 1200.0);
    port.setGeometry(0.1, 0.03, 1.0, 0.0);

    ExhaustNetwork net;
    ExhaustNetwork::PrimarySpec spec;
    spec.port = &port;
    spec.length = 0.5;
    spec.diameter = 0.04;
    spec.orificeArea = 5e-4;
    net.addPrimary(spec);
    Pipe1D::Section s;
    s.length = 0.3;
    s.inletDiameter = s.outletDiameter = 0.04;
    ExhaustNetwork::Parameters np;
    np.collectorSections.push_back(s);
    np.cellSize = 0.02;
    np.wallTemperature = Tatm;   // no wall heat exchange
    net.build(np);

    const double e0 = port.kineticEnergy() + net.totalEnergy();
    const double m0 = port.mass() + net.totalMass();
    const double dt = 1e-5;
    for (int i = 0; i < 200; ++i) net.exchangePort(0, dt);   // 2 ms of blowdown, pipe frozen
    const double e1 = port.kineticEnergy() + net.totalEnergy();
    const double m1 = port.mass() + net.totalMass();
    EXPECT_LT(port.pressure(), 2.9 * Patm);     // something actually flowed
    EXPECT_NEAR(m1, m0, 1e-9 * m0);
    EXPECT_NEAR(e1, e0, 1e-6 * e0);
}
