// Network-level tests for the 1-D exhaust model: conservation through the
// port coupling and the whole network.
//
// Copyright (c) 2026 Arcindweld. MIT License (see LICENSE).

#include <gtest/gtest.h>

#include "../include/exhaust_network.h"
#include "../include/gas_system.h"
#include "../include/units.h"

#include <algorithm>
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

    // returns the mass (kg) added to the ports by this call. The feed is a fixed molar rate, or with `smooth` a
    // sin^2 pulse of the same mean (a valve-like opening and closing)
    double drive(double t, double dt, bool smooth = false) {
        const double period = 2.0 / (6000.0 / 60.0);   // one engine cycle, s
        double added = 0.0;
        for (int i = 0; i < 4; ++i) {
            const double phase = std::fmod(t + i * period / 4, period) / period;
            if (phase < 0.12) {
                // feed hot gas (like a valve opening)
                const double sn = std::sin(3.14159265358979323846 * phase / 0.12);
                const double nDot = smooth ? 1.0 * sn * sn : 0.5;   // mol/s
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

// Give a port bulk motion along +x the way the engine does: gas pushed in from a reservoir through
// GasSystem::flow, which hands the receiving system the momentum of the incoming fraction.
void spinUp(GasSystem &port, double reservoirPressure, double reservoirTemperature) {
    GasSystem reservoir;
    reservoir.initialize(reservoirPressure, 10e-3, reservoirTemperature);
    GasSystem::FlowParameters fp;
    fp.k_flow = GasSystem::k_28inH2O(150.0);
    fp.dt = 1e-5;
    fp.direction_x = 1.0;
    fp.direction_y = 0.0;
    fp.crossSectionArea_0 = 0.0;
    fp.crossSectionArea_1 = 0.1 * 0.03;
    fp.system_0 = &reservoir;
    fp.system_1 = &port;
    for (int i = 0; i < 20; ++i) GasSystem::flow(fp);
}

// One port emptying into (or filling from) a 0.5 m primary and a 2 m collector, no wall friction or heat exchange,
// for 1.5 ms (the waves stay clear of the open end); returns the relative change of (port total energy + pipe
// energy) and of the mass.
struct ExchangeBalance { double energy, mass, portVelocity, portPressure, massOut; };
ExchangeBalance portBlowdown(GasSystem &port) {
    ExhaustNetwork net;
    ExhaustNetwork::PrimarySpec spec;
    spec.port = &port;
    spec.length = 0.5;
    spec.diameter = 0.04;
    spec.orificeArea = 5e-4;
    net.addPrimary(spec);
    Pipe1D::Section s;
    s.length = 2.0;
    s.inletDiameter = s.outletDiameter = 0.04;
    ExhaustNetwork::Parameters np;
    np.collectorSections.push_back(s);
    np.cellSize = 0.01;
    np.frictionFactor = 0.0;   // no wall friction, and with it no wall heat exchange (Reynolds analogy)
    net.build(np);

    ExchangeBalance res;
    res.portVelocity = port.velocity_x();
    const double e0 = port.totalEnergy() + net.totalEnergy();
    const double m0 = port.mass() + net.totalMass();
    for (int i = 0; i < 150; ++i) net.advance(1e-5);
    net.finishStep();
    res.energy = (port.totalEnergy() + net.totalEnergy() - e0) / e0;
    res.mass = (port.mass() + net.totalMass() - m0) / m0;
    res.portPressure = port.pressure();
    res.massOut = net.massOut();
    return res;
}

// The pulsed four-port network driven the way the engine drives it: per fluid sub-step the ports are fed, then
// the network advances (its port boundaries move gas between the ports and the pipes); once per simulation step
// the step is closed. Returns the mean and RMS port gauge pressure and the forward flow per cycle over the last
// 10 of 20 cycles at 6000 rpm, the share of it that came back, and how often the port flow changed direction.
struct CoupledRun { double portMean, portRms, forwardPerCycle, backflowShare, reversalsPerCycle, pipeSteps; };
CoupledRun runCoupled(double cellSize, double simulationFrequency, bool holdPipes = true) {
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
    ExhaustNetwork::Parameters np = fourIntoOne();
    np.cellSize = cellSize;
    np.holdPipes = holdPipes;
    net.build(np);

    const int fluidSteps = 8;
    const double dt = 1.0 / simulationFrequency, dtFluid = dt / fluidSteps;
    const double period = 0.02, tAverage = 0.2, tEnd = 0.4;
    double t = 0, pSum = 0, p2Sum = 0, tSum = 0, forward = 0, back = 0, lastSign = 0;
    long reversals = 0;
    double pipeSteps = 0, stepCount = 0;
    while (t < tEnd - 1e-12) {
        for (int s = 0; s < fluidSteps; ++s) {
            ports.drive(t + s * dtFluid, dtFluid, true);
            const double m0 = ports.port[0].mass();
            net.advance(dtFluid);
            if (t >= tAverage) {
                // port 0's exchange over this fluid sub-step (the network may have sub-stepped)
                const double m = m0 - ports.port[0].mass();
                if (m > 0) forward += m; else back -= m;
                if (m != 0) {
                    const double sign = m > 0 ? 1 : -1;
                    if (lastSign != 0 && sign != lastSign) ++reversals;
                    lastSign = sign;
                }
                const double gp = ports.port[0].pressure() - Patm;
                pSum += gp * dtFluid;
                p2Sum += gp * gp * dtFluid;
                tSum += dtFluid;
            }
        }
        net.finishStep();
        if (t >= tAverage) { pipeSteps += net.lastSubSteps(); ++stepCount; }
        t += dt;
    }
    const double cycles = (tEnd - tAverage) / period;
    CoupledRun r;
    r.portMean = pSum / tSum;
    r.portRms = std::sqrt(std::max(0.0, p2Sum / tSum - r.portMean * r.portMean));
    r.forwardPerCycle = forward / cycles;
    r.backflowShare = back / std::max(forward, 1e-30);
    r.reversalsPerCycle = reversals / cycles;
    r.pipeSteps = pipeSteps / std::max(stepCount, 1.0);
    return r;
}

} // namespace

// (d) Mass conservation over 1 s of pulsed running: mass added to the ports
// = mass left in the ports + mass in the pipes + mass out of the open end, to
// round-off (1e-9 of the throughput asserted).
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

    // driven like the engine: 10 kHz steps of 8 fluid sub-steps; the ports are fed, then the network advances
    const double dt = 1e-4, dtFluid = dt / 8;
    double added = 0.0, t = 0.0;
    int steps = 0;
    while (t < 1.0) {
        for (int s = 0; s < 8; ++s) {
            added += ports.drive(t + s * dtFluid, dtFluid);
            net.advance(dtFluid);
        }
        net.finishStep();
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
    EXPECT_NEAR(balance / added, 0.0, 1e-9);
    EXPECT_NEAR(net.massIn(), added - (portMass - portMass0), 1e-9 * added);
}

// Energy conservation of the port boundary. Total energy (port internal + bulk kinetic energy, plus the pipes'
// total energy) is conserved: gas leaving the port takes its enthalpy and its share of the port's bulk kinetic
// energy into the pipe as total energy; gas coming back from the pipe brings its stagnation enthalpy, and the
// port's bulk motion that it slows turns into heat. Four cases: a hot port at rest emptying into the pipe, the
// same port with bulk motion along the runner (as the engine's exhaust runner has), and a port below atmospheric
// pressure filling from the pipe, at rest and moving. (The first version kept the leaving gas's momentum in the
// port, which made kinetic energy out of nothing, and dropped the bulk kinetic energy on back-flow.)
TEST(ExhaustNetworkTests, PortExchangeConservesEnergy) {
    for (int moving = 0; moving < 2; ++moving) {
        for (int filling = 0; filling < 2; ++filling) {
            GasSystem port;
            const double p = filling ? 0.5 * Patm : 3.0 * Patm, T = filling ? 400.0 : 1200.0;
            port.initialize(p, 100e-6, T);
            port.setGeometry(0.1, 0.03, 1.0, 0.0);
            if (moving) spinUp(port, 2.0 * p, T);
            const double p0 = port.pressure();
            const ExchangeBalance b = portBlowdown(port);
            std::cout << "[  values  ] " << (filling ? "filling" : "emptying") << ", port at "
                      << b.portVelocity << " m/s: energy change " << b.energy << ", mass change " << b.mass
                      << ", port " << p0 << " -> " << b.portPressure << " Pa\n";
            if (moving) EXPECT_GT(std::fabs(b.portVelocity), 10.0) << "the spin-up gave the port no motion";
            if (filling) EXPECT_GT(b.portPressure, 1.2 * p0);   // something actually flowed
            else EXPECT_LT(b.portPressure, 0.8 * p0);
            EXPECT_LT(std::fabs(b.massOut), 1e-15);             // nothing reached the open end
            EXPECT_NEAR(b.mass, 0.0, 1e-12);
            EXPECT_NEAR(b.energy, 0.0, 1e-12);
        }
    }
}

// Steady flow through the port boundary: a port held at 1.3 bar, 900 K (a 1 m^3 volume) feeding a frictionless
// 4-into-1 whose collector is open: after the start-up waves have died away the pipes sit at atmospheric
// pressure and the mass flow is the orifice's from the port state to 1 atm, the face velocity is that flow
// carried by port gas at 1 atm.
TEST(ExhaustNetworkTests, SteadyPortFlowIsTheOrificeFlow) {
    GasSystem port;
    port.initialize(1.3e5, 1.0, 900.0);
    ExhaustNetwork net;
    ExhaustNetwork::PrimarySpec spec;
    spec.port = &port;
    spec.length = 0.4;
    spec.diameter = 0.04;
    spec.orificeArea = 4e-4;
    net.addPrimary(spec);
    ExhaustNetwork::Parameters np;
    Pipe1D::Section s;
    s.length = 0.4;
    s.inletDiameter = s.outletDiameter = 0.04;
    np.collectorSections.push_back(s);
    np.cellSize = 0.01;
    np.frictionFactor = 0.0;
    net.build(np);
    for (int i = 0; i < 6000; ++i) {      // 0.3 s
        net.advance(5e-5);
        net.finishStep();
        port.reset(1.3e5, 900.0, port.mix());
    }
    const double expected = net.orificeFlow(4e-4, 1.3e5, 900.0, Patm);
    std::cout << "[  values  ] port flow " << net.portMassFlow(0) << " kg/s, outlet " << net.outletMassFlow()
              << " kg/s, orifice to 1 atm " << expected << " kg/s\n";
    EXPECT_NEAR(net.portMassFlow(0), expected, 0.005 * expected);
    EXPECT_NEAR(net.outletMassFlow(), expected, 0.005 * expected);
}

// The coupled port / pipe solution converges and does not chatter. The pulsed four-port network with a
// valve-like (sin^2) feed, driven as the engine drives it: the port flow changes direction a few times per cycle
// (the waves), not on every sub-step; the port pressure history (RMS of the gauge pressure), the forward flow and
// the back-flow hardly depend on the simulation frequency (0.2 % between 10 and 40 kHz on 1 cm cells), and approach
// the fine-grid values as the cells shrink (RMS -23 % at 4 cm, -8 % at 2 cm against 1 cm; 1 cm is within 1 % of
// 0.5 cm).
// The port's own volume against the pipe rings near c A / (2 pi V), about 2 kHz for this 60 cc port: above what 4 cm
// cells carry (docs/exhaust-1d.md, section 5), hence the coarse-grid error. With the orifice computed from the
// first cell's state and added to the cell (the first coupling) the flow reversed on 90 % of the sub-steps, the
// forward flow was 6-19 times the gas fed, and both moved with the grid and with the frequency without settling.
TEST(ExhaustNetworkTests, CoupledPortsConvergeWithCellSizeAndFrequency) {
    const double cells[3] = { 0.04, 0.02, 0.01 };
    CoupledRun run[3];
    for (int k = 0; k < 3; ++k) run[k] = runCoupled(cells[k], 20000);
    const CoupledRun f10 = runCoupled(0.01, 10000), f40 = runCoupled(0.01, 40000);
    auto show = [](const char *tag, const CoupledRun &r) {
        std::cout << "[  values  ] " << tag << ": port mean " << r.portMean << " Pa, rms " << r.portRms
                  << " Pa, forward " << 1e6 * r.forwardPerCycle << " mg/cycle, back-flow " << 100 * r.backflowShare
                  << " %, reversals " << r.reversalsPerCycle << " /cycle\n";
    };
    show("4 cm 20 kHz", run[0]); show("2 cm 20 kHz", run[1]); show("1 cm 20 kHz", run[2]);
    show("1 cm 10 kHz", f10); show("1 cm 40 kHz", f40);
    const CoupledRun &ref = run[2];
    ASSERT_GT(ref.portRms, 1000.0);
    // no chatter: 1600 (10 kHz) to 6400 (40 kHz) sub-steps per cycle
    const CoupledRun *all[] = { &run[0], &run[1], &run[2], &f10, &f40 };
    for (const CoupledRun *r : all) EXPECT_LT(r->reversalsPerCycle, 60.0);
    // the time step: 10 and 40 kHz against 20 kHz on 1 cm cells
    for (const CoupledRun *r : { &f10, &f40 }) {
        EXPECT_NEAR(r->portRms, ref.portRms, 0.01 * ref.portRms);
        EXPECT_NEAR(r->forwardPerCycle, ref.forwardPerCycle, 0.01 * ref.forwardPerCycle);
        EXPECT_NEAR(r->backflowShare, ref.backflowShare, 0.01);
    }
    // the grid: each halving of the cells at least halves the distance to the 1 cm result
    auto converging = [&](double coarse, double mid, double fine) {
        return std::fabs(mid - fine) < 0.5 * std::fabs(coarse - fine);
    };
    EXPECT_TRUE(converging(run[0].portRms, run[1].portRms, ref.portRms));
    EXPECT_TRUE(converging(run[0].forwardPerCycle, run[1].forwardPerCycle, ref.forwardPerCycle));
    EXPECT_TRUE(converging(run[0].backflowShare, run[1].backflowShare, ref.backflowShare));
    EXPECT_NEAR(run[1].portRms, ref.portRms, 0.15 * ref.portRms);
}

// Multi-rate stepping (Parameters::holdPipes): the pipes at their own CFL step while the ports keep the fluid
// sub-step against the pipes' held state, against lock-step (the pipes advanced with every sub-step) on the coupled
// ports at 20 kHz. On 4 and 2 cm cells lock-step runs the pipes at a small fraction of their CFL step, where the
// upwind scheme is most diffusive; held, they run near their CFL step and come out closer to the 1 cm result. On
// 1 cm cells at 40 kHz, where both resolve the waves, the two agree.
TEST(ExhaustNetworkTests, HeldPipesConvergeToTheFineGrid) {
    const CoupledRun ref = runCoupled(0.01, 20000, false);   // CFL-limited: lock-step either way
    for (double cell : { 0.04, 0.02 }) {
        const CoupledRun held = runCoupled(cell, 20000, true), lock = runCoupled(cell, 20000, false);
        std::cout << "[  values  ] " << 100 * cell << " cm 20 kHz, held / lock-step / 1 cm: rms " << held.portRms
                  << " / " << lock.portRms << " / " << ref.portRms << " Pa, forward " << 1e6 * held.forwardPerCycle
                  << " / " << 1e6 * lock.forwardPerCycle << " / " << 1e6 * ref.forwardPerCycle
                  << " mg/cycle, back-flow " << 100 * held.backflowShare << " / " << 100 * lock.backflowShare << " / "
                  << 100 * ref.backflowShare << " %, pipe steps per step " << held.pipeSteps << " / "
                  << lock.pipeSteps << "\n";
        EXPECT_LT(held.pipeSteps, 0.5 * lock.pipeSteps);   // the point of it
        EXPECT_LT(std::fabs(held.portRms - ref.portRms), std::fabs(lock.portRms - ref.portRms));
        EXPECT_LT(std::fabs(held.forwardPerCycle - ref.forwardPerCycle),
                  std::fabs(lock.forwardPerCycle - ref.forwardPerCycle));
        EXPECT_LT(std::fabs(held.backflowShare - ref.backflowShare),
                  std::fabs(lock.backflowShare - ref.backflowShare));
    }
    const CoupledRun held = runCoupled(0.01, 40000, true), lock = runCoupled(0.01, 40000, false);
    std::cout << "[  values  ] 1 cm 40 kHz, held / lock-step: rms " << held.portRms << " / " << lock.portRms
              << " Pa, forward " << 1e6 * held.forwardPerCycle << " / " << 1e6 * lock.forwardPerCycle
              << " mg/cycle, pipe steps per step " << held.pipeSteps << " / " << lock.pipeSteps << "\n";
    EXPECT_LT(held.pipeSteps, 0.5 * lock.pipeSteps);
    EXPECT_NEAR(held.portRms, lock.portRms, 0.01 * lock.portRms);
    EXPECT_NEAR(held.forwardPerCycle, lock.forwardPerCycle, 0.01 * lock.forwardPerCycle);
    EXPECT_NEAR(held.backflowShare, lock.backflowShare, 0.01);
}
