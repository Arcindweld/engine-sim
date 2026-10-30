// Network of 1-D pipes for one exhaust system: N primaries that join in a
// constant-pressure junction, then a collector run (cones and steps allowed)
// that ends open to the atmosphere.
//
// Copyright (c) 2026 Arcindweld. MIT License (see LICENSE).
//
// Coupling to the rest of engine-sim: each primary's inlet is fed from the
// cylinder's exhaust-port GasSystem (the lumped runner volume that the cylinder
// already flows into) through a compressible orifice whose effective area is
// derived from the script's `primary_flow_rate`.  The port exchange runs at the
// fluid sub-step rate (it only needs the pipe's end state); the pipes advance
// once per simulation step with their own CFL sub-stepping, driven by the
// port fluxes averaged over that step.  Mass and energy that leave a port
// GasSystem arrive in the pipe exactly.

#ifndef ATG_ENGINE_SIM_EXHAUST_NETWORK_H
#define ATG_ENGINE_SIM_EXHAUST_NETWORK_H

#include "pipe_1d.h"
#include "gas_system.h"

#include <vector>

class ExhaustNetwork {
    public:
        struct Parameters {
            std::vector<Pipe1D::Section> collectorSections;  // from the junction to the open end
            double cellSize = 0.02;
            int maxCellsPerPipe = 200;
            double wallTemperature = 700.0;
            double frictionFactor = 0.005;
            double atmosphericPressure = 101325.0;
            double atmosphericTemperature = 298.15;
            double initialTemperature = 298.15;
            double cfl = 0.8;
            double molarMass = 0.02897;         // kg/mol of the GasSystem gas
            int degreesOfFreedom = 5;           // GasSystem gas model -> gamma = 1 + 2/dof
        };

        struct PrimarySpec {
            GasSystem *port = nullptr;   // lumped exhaust-port volume of the cylinder
            double length = 0.0;         // m, port exit to the junction
            double diameter = 0.0;       // m
            double orificeArea = 0.0;    // m^2, effective area Cd*A of the port -> pipe entry
        };

        struct Primary {
            PrimarySpec spec;
            Pipe1D pipe;
            double accumulated[Pipe1D::FluxCount];   // advective flux x time since the last step
            double averaged[Pipe1D::FluxCount];      // the same per unit time, used during the step
            double jetSpeed;                          // largest jet speed seen since the last step
            double lastMassFlow;                      // kg/s into the pipe, last exchange
        };

        struct Junction {
            double pressure = 0.0;       // last solution, used as the Newton start
        };

    public:
        ExhaustNetwork();
        ~ExhaustNetwork();

        int addPrimary(const PrimarySpec &spec);
        void build(const Parameters &params);
        inline bool isBuilt() const { return m_built; }

        // Per fluid sub-step: move gas between the port GasSystem and the pipe end.
        void exchangePort(int primary, double dt);
        // Once per simulation step.
        void step(double dt);

        // Outlet (open end) quantities, averaged over the last step
        inline double outletMassFlow() const { return m_outletMassFlow; }                 // kg/s, positive out
        inline double outletMassFlowDerivative() const { return m_outletMassFlowDerivative; } // kg/s^2
        // Far-field monopole pressure at 1 m: rho0/(4 pi r) dQ/dt with Q = mdot/rho0
        inline double radiatedPressure() const { return m_outletMassFlowDerivative / (4 * 3.14159265358979323846); }
        double outletPressure() const;
        double outletTemperature() const;
        GasSystem::Mix outletMix() const;

        // Accounting
        double totalMass() const;
        double totalEnergy() const;
        inline double massIn() const { return m_massIn; }       // through the ports, net into the pipes
        inline double massOut() const { return m_massOut; }     // through the open end, net out
        inline double initialMass() const { return m_initialMass; }
        inline int lastSubSteps() const { return m_lastSubSteps; }
        int cellCount() const;

        inline int primaryCount() const { return static_cast<int>(m_primaries.size()); }
        inline Pipe1D &primaryPipe(int i) { return m_primaries[i].pipe; }
        inline const Pipe1D &primaryPipe(int i) const { return m_primaries[i].pipe; }
        inline int collectorPipeCount() const { return static_cast<int>(m_collector.size()); }
        inline Pipe1D &collectorPipe(int i) { return *m_collector[i]; }
        inline const Pipe1D &collectorPipe(int i) const { return *m_collector[i]; }
        inline const Parameters &parameters() const { return m_params; }
        inline double junctionPressure() const { return m_junction.pressure; }

        // Compressible orifice mass flow (kg/s) from stagnation (p0, T0) to pd
        double orificeFlow(double CdA, double p0, double T0, double pd) const;
        double jetSpeed(double p0, double T0, double pd) const;

        // Characteristic junctions (exposed for tests)
        // Two pipes meeting end to end: a's right end to b's left end.
        void junction2(Pipe1D &a, Pipe1D &b, Junction &j) const;
        struct JunctionEnd {
            Pipe1D *pipe;
            int side;
        };
        void junctionN(JunctionEnd *ends, int count, Junction &j) const;

    private:
        void applyBoundaries();
        void applyPortFlux(Primary &prim);
        double gammaValue() const { return m_gamma; }

    private:
        Parameters m_params;
        std::vector<Primary> m_primaries;
        std::vector<Pipe1D *> m_collector;
        std::vector<Junction> m_collectorJunctions;
        Junction m_junction;
        std::vector<JunctionEnd> m_junctionEnds;
        bool m_built;

        double m_gamma, m_R, m_cp;
        double m_atmScalars[Pipe1D::ScalarCount];

        double m_outletMassFlow;
        double m_outletMassFlowDerivative;
        double m_massIn, m_massOut, m_initialMass;
        int m_lastSubSteps;
};

#endif /* ATG_ENGINE_SIM_EXHAUST_NETWORK_H */
