// One-dimensional compressible duct solver for exhaust pipes.
//
// Copyright (c) 2026 Arcindweld. MIT License (see LICENSE).
//
// Quasi-1-D Euler equations for an ideal gas in a duct of slowly varying
// cross-section, with wall friction (Fanning) and Reynolds-analogy wall heat
// transfer, solved with a MUSCL-Hancock finite-volume scheme and an HLLC
// Riemann solver (Toro, "Riemann Solvers and Numerical Methods for Fluid
// Dynamics", 3rd ed., ch. 10 and 14).  The scheme, the open-end ghost state
// and the characteristic junctions are a port of the validated JavaScript
// reference `wsFactory()` in the ExhaustTuner calculator (same author).
//
// Conserved variables per cell (all multiplied by the cell area A):
//   q1 = rho A, q2 = rho u A, q3 = E A with E = p/(gamma-1) + rho u^2 / 2,
//   qs[k] = rho A x_k for ScalarCount passive scalars (mixture fractions).
//
// The pipe does not step itself: the owning network calls computePrimitives(),
// fills the two end fluxes (open end, wall, fed port or junction) and then
// calls advance(dt) with a dt that satisfies the CFL condition.

#ifndef ATG_ENGINE_SIM_PIPE_1D_H
#define ATG_ENGINE_SIM_PIPE_1D_H

#include <cmath>
#include <vector>

class Pipe1D {
    public:
        static constexpr int ScalarCount = 2;      // x[0] = fuel fraction, x[1] = O2 fraction
        static constexpr int FluxCount = 3 + ScalarCount;

        struct Section {
            double length = 0.0;           // m
            double inletDiameter = 0.0;    // m, at the start of the section
            double outletDiameter = 0.0;   // m, at the end (linear taper in between)
        };

        struct Parameters {
            std::vector<Section> sections;  // contiguous; a diameter jump belongs between two pipes
            double cellSize = 0.04;         // m, target; the actual size is length / cellCount
            int minCells = 6;
            int maxCells = 200;             // long pipes are coarsened so a run never exceeds this
            double wallTemperature = 700.0; // K
            double initialPressure = 101325.0;
            double initialTemperature = 298.15;
            double initialScalars[ScalarCount] = { 0.0, 0.0 };
            double frictionFactor = 0.005;  // Fanning friction factor f (tau_w = f rho u^2 / 2)
            double heatCapacityRatio = 1.4;
            double specificGasConstant = 287.0; // J/(kg K)
        };

        // Riemann state on one side of a face
        struct State {
            double rho, u, p;
            double x[ScalarCount];
        };

    public:
        Pipe1D();
        ~Pipe1D();

        void initialize(const Parameters &params);

        // Riemann solver and gas relations (shared with the junctions)
        void hllc(const State &L, const State &R, double *flux) const;
        inline double gamma() const { return m_gamma; }
        inline double gasConstant() const { return m_R; }
        inline double cp() const { return m_cp; }
        inline double soundSpeed(double p, double rho) const;

        // Geometry
        inline int cellCount() const { return m_n; }
        inline double cellSize() const { return m_dx; }
        inline double length() const { return m_n * m_dx; }
        inline double cellArea(int i) const { return m_A[i]; }
        inline double faceArea(int j) const { return m_Af[j]; }   // j in [0, n]
        inline double endFaceArea(int side) const { return side > 0 ? m_Af[m_n] : m_Af[0]; }
        inline int endCell(int side) const { return side > 0 ? m_n - 1 : 0; }
        inline int endFace(int side) const { return side > 0 ? m_n : 0; }
        inline double wallTemperature() const { return m_Tw; }
        inline void setWallTemperature(double T) { m_Tw = T; }

        // Primitive state (valid after computePrimitives())
        void computePrimitives();
        inline double density(int i) const { return m_r[i]; }
        inline double velocity(int i) const { return m_u[i]; }
        inline double pressure(int i) const { return m_p[i]; }
        inline double temperature(int i) const { return m_p[i] / (m_R * m_r[i]); }
        inline double scalar(int i, int k) const { return m_x[k][i]; }
        State state(int i) const;
        // max |u| + c over the cells, cached by computePrimitives() and advance()
        inline double maxSignalSpeed() const { return m_maxSignalSpeed; }

        // End fluxes (per face, already multiplied by the face area): mass,
        // momentum, energy, scalars.  side = -1 is the left end (x = 0),
        // side = +1 the right end (x = L).
        void setEndFlux(int side, const double *flux);
        void getEndFlux(int side, double *flux) const;
        void applyWallEnd(int side);
        // Open end to still gas at (pAtm, TAtm) with the given scalar fractions
        void applyOpenEnd(int side, double pAtm, double TAtm, const double *scalarsAtm);

        // One MUSCL-Hancock step.  Requires the primitive state (from
        // initialize(), computePrimitives() or the previous advance()) and both
        // end fluxes for the current state.  Leaves the primitive state and
        // maxSignalSpeed() updated for the new state.
        void advance(double dt);

        // Diagnostics
        double totalMass() const;
        double totalEnergy() const;
        // Mass flow through an end face in the +x direction, from the last advance() (kg/s)
        inline double endMassFlow(int side) const { return m_f[0][endFace(side)]; }
        // Mass flow, time-averaged over the last network step (set by the owner)
        inline double averagedEndMassFlow(int side) const { return side > 0 ? m_avgMassFlow[1] : m_avgMassFlow[0]; }
        inline void setAveragedEndMassFlow(int side, double value) { (side > 0 ? m_avgMassFlow[1] : m_avgMassFlow[0]) = value; }

        // Direct state manipulation for tests and initial conditions
        void setCellState(int i, double rho, double u, double p, const double *scalars);
        // Add mass (kg), momentum (kg m/s), energy (J) and scalar mass (kg) to
        // one cell and refresh its primitive state (used by the port coupling).
        void addToCell(int i, double dMass, double dMomentum, double dEnergy, const double *dScalarMass);
        inline double cellMass(int i) const { return m_q1[i] * m_dx; }

    private:
        static inline double minmod(double a, double b) {
            if (a * b <= 0) return 0.0;
            return (std::fabs(a) < std::fabs(b)) ? a : b;
        }

    private:
        int m_n;
        double m_dx;
        double m_gamma, m_R, m_cv, m_cp;
        double m_Tw;
        double m_friction;

        std::vector<double> m_A, m_Af, m_D;
        std::vector<double> m_invA, m_invD, m_invR;
        double m_invDx, m_invG1;
        std::vector<double> m_q1, m_q2, m_q3;
        std::vector<double> m_qs[ScalarCount];
        std::vector<double> m_r, m_u, m_p;
        std::vector<double> m_x[ScalarCount];
        std::vector<State> m_bL, m_bR;           // evolved face states of each cell
        std::vector<double> m_f[FluxCount];      // face fluxes, size n + 1
        double m_maxSignalSpeed;
        double m_avgMassFlow[2];
};

inline double Pipe1D::soundSpeed(double p, double rho) const {
    return std::sqrt(m_gamma * p / rho);
}

#endif /* ATG_ENGINE_SIM_PIPE_1D_H */
