#include "../include/exhaust_system.h"

#include "../include/units.h"

#include <algorithm>
#include <cmath>

ExhaustSystem::ExhaustSystem() {
    m_primaryFlowRate = 0;
    m_outletFlowRate = 0;
    m_collectorCrossSectionArea = 0;
    m_length = 0;
    m_primaryTubeLength = 0;
    m_audioVolume = 0;
    m_velocityDecay = 0;
    m_flow = 0;
    m_index = -1;
    m_impulseResponse = nullptr;
    m_legacy = false;
}

ExhaustSystem::~ExhaustSystem() {
    /* void */
}

void ExhaustSystem::initialize(const Parameters &params) {
    const double systemWidth = std::sqrt(params.collectorCrossSectionArea);
    const double volume = params.collectorCrossSectionArea * params.length;
    const double systemLength = params.length;
    m_system.initialize(
            units::pressure(1.0, units::atm),
            volume,
            units::celcius(25.0));
    m_system.setGeometry(
        systemLength,
        systemWidth,
        1.0,
        0.0);

    m_atmosphere.initialize(
        units::pressure(1.0, units::atm),
        units::volume(1000.0, units::m3),
        units::celcius(25.0));
    m_atmosphere.setGeometry(
        units::distance(10.0, units::m),
        units::distance(10.0, units::m),
        1.0,
        0.0);

    m_primaryFlowRate = params.primaryFlowRate;
    m_audioVolume = params.audioVolume;
    m_outletFlowRate = params.outletFlowRate;
    m_collectorCrossSectionArea = params.collectorCrossSectionArea;
    m_velocityDecay = params.velocityDecay;
    m_impulseResponse = params.impulseResponse;
    m_length = params.length;
    m_primaryTubeLength = params.primaryTubeLength;

    m_parameters = params;
    m_legacy = params.selectsLegacy();
}

bool ExhaustSystem::Parameters::hasGeometry() const {
    return primaryDiameter > 0
        || collectorDiameter > 0
        || collectorLength >= 0
        || megaphoneLength > 0
        || megaphoneExitDiameter > 0
        || tailpipeLength > 0
        || tailpipeDiameter > 0
        || pipeCellSize > 0;
}

bool ExhaustSystem::Parameters::selectsLegacy() const {
    if (legacy) return true;
    if (exhaustModel == "legacy") return true;
    if (exhaustModel == "1d") return false;
    // "auto" (and anything unrecognised): the 1-D model only when the file
    // describes 1-D geometry
    return !hasGeometry();
}

void ExhaustSystem::destroy() {
    /* void */
}

void ExhaustSystem::process(double dt) {
    if (!m_legacy) return;

    GasSystem::Mix airMix;
    airMix.p_fuel = 0;
    airMix.p_inert = 1.0;
    airMix.p_o2 = 0.0;

    m_atmosphere.reset(units::pressure(1.0, units::atm), units::celcius(25.0), airMix);
    GasSystem::FlowParameters flowParams;
    flowParams.crossSectionArea_0 = m_collectorCrossSectionArea;
    flowParams.crossSectionArea_1 = units::area(10, units::m2);
    flowParams.direction_x = 1.0;
    flowParams.direction_y = 0.0;
    flowParams.dt = dt;
    flowParams.system_0 = &m_atmosphere;
    flowParams.system_1 = &m_system;
    flowParams.k_flow = m_outletFlowRate;

    m_flow = m_system.flow(flowParams);

    m_system.dissipateExcessVelocity();
    m_system.updateVelocity(dt, m_velocityDecay);
}

int ExhaustSystem::attachPrimary(GasSystem *port, double length, double portArea, double kFlow) {
    const double pi = 3.14159265358979323846;
    ExhaustNetwork::PrimarySpec spec;
    spec.port = port;
    spec.length = length;
    spec.diameter = (m_parameters.primaryDiameter > 0)
        ? m_parameters.primaryDiameter
        : std::sqrt(4 * portArea / pi);
    const double pipeArea = pi * spec.diameter * spec.diameter / 4;
    // GasSystem::flowRate gives n_dot = k * p0 sqrt(gamma / (R T0)) * f(pr); the
    // physical orifice gives m_dot = Cd A p0 sqrt(gamma M / (R T0)) f(pr), so
    // Cd A = k sqrt(M).  Never wider than the port or the pipe it feeds.
    const double orificeArea = kFlow * std::sqrt(units::AirMolecularMass);
    spec.orificeArea = std::min(orificeArea, std::min(portArea, pipeArea));
    return m_network.addPrimary(spec);
}

void ExhaustSystem::buildNetwork() {
    if (m_legacy || m_network.isBuilt()) return;

    const double pi = 3.14159265358979323846;
    const Parameters &p = m_parameters;

    ExhaustNetwork::Parameters np;
    np.cellSize = (p.pipeCellSize > 0) ? p.pipeCellSize : DefaultPipeCellSize;
    np.wallTemperature = p.wallTemperature;
    np.atmosphericPressure = units::pressure(1.0, units::atm);
    np.atmosphericTemperature = units::celcius(25.0);
    np.initialTemperature = units::celcius(25.0);
    np.molarMass = units::AirMolecularMass;
    np.degreesOfFreedom = m_system.degreesOfFreedom();

    const double collectorDiameter = (p.collectorDiameter > 0)
        ? p.collectorDiameter
        : std::sqrt(4 * p.collectorCrossSectionArea / pi);
    const double collectorLength = (p.collectorLength >= 0) ? p.collectorLength : p.length;

    Pipe1D::Section s;
    double lastDiameter = collectorDiameter;
    if (collectorLength > 0) {
        s.length = collectorLength;
        s.inletDiameter = s.outletDiameter = collectorDiameter;
        np.collectorSections.push_back(s);
    }
    if (p.megaphoneLength > 0 && p.megaphoneExitDiameter > 0) {
        s.length = p.megaphoneLength;
        s.inletDiameter = lastDiameter;
        s.outletDiameter = p.megaphoneExitDiameter;
        np.collectorSections.push_back(s);
        lastDiameter = p.megaphoneExitDiameter;
    }
    if (p.tailpipeLength > 0) {
        s.length = p.tailpipeLength;
        s.inletDiameter = s.outletDiameter = (p.tailpipeDiameter > 0) ? p.tailpipeDiameter : lastDiameter;
        np.collectorSections.push_back(s);
        lastDiameter = s.outletDiameter;
    }

    m_network.build(np);
}

void ExhaustSystem::exchangePort(int primaryIndex, double dt) {
    if (!m_network.isBuilt()) buildNetwork();
    m_network.exchangePort(primaryIndex, dt);
}

void ExhaustSystem::processNetwork(double dt) {
    if (m_legacy) return;
    if (!m_network.isBuilt()) buildNetwork();
    m_network.step(dt);

    // Mirror the outlet state into the legacy GasSystem so that readers of
    // getSystem() (exhaust O2 gauge, average output signal) keep working.
    m_system.reset(
        m_network.outletPressure(),
        m_network.outletTemperature(),
        m_network.outletMix());
    m_flow = m_network.outletMassFlow() * dt / units::AirMolecularMass;
}
