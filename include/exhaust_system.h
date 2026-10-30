#ifndef ATG_ENGINE_SIM_EXHAUST_SYSTEM_H
#define ATG_ENGINE_SIM_EXHAUST_SYSTEM_H

#include "part.h"

#include "gas_system.h"
#include "impulse_response.h"
#include "exhaust_network.h"

#include <string>
#include <vector>

class ExhaustSystem : public Part {
    friend class Engine;

    public:
        struct Parameters {
            double length;
            double collectorCrossSectionArea;
            double outletFlowRate;
            double primaryTubeLength;
            double primaryFlowRate;
            double velocityDecay;
            double audioVolume;
            ImpulseResponse *impulseResponse;

            // 1-D gas-dynamics model (see docs/exhaust-1d.md).  All optional;
            // zero or negative means "derive from the legacy parameters".
            // exhaustModel: "auto" (1-D only when one of the geometry inputs
            // below is set, legacy otherwise), "1d" or "legacy".
            std::string exhaustModel = "auto";
            bool legacy = false;              // alias for exhaustModel = "legacy"
            double primaryDiameter = 0.0;     // default: head exhaust_runner_cross_section_area
            double collectorDiameter = 0.0;   // default: collector_cross_section_area
            double collectorLength = -1.0;    // default: length
            double megaphoneLength = 0.0;     // cone after the collector (0 = none)
            double megaphoneExitDiameter = 0.0;
            double tailpipeLength = 0.0;      // straight pipe after the cone (0 = none)
            double tailpipeDiameter = 0.0;    // default: the diameter it follows
            double pipeCellSize = 0.0;        // m; 0 = default (DefaultPipeCellSize)
            double wallTemperature = 700.0;   // K

            // true when any of the 1-D geometry inputs differs from its default
            bool hasGeometry() const;
            // the model that these parameters select
            bool selectsLegacy() const;
        };

        static constexpr double DefaultPipeCellSize = 0.04;   // m (see docs/exhaust-1d.md, performance)

    public:
        ExhaustSystem();
        virtual ~ExhaustSystem();

        void initialize(const Parameters &params);
        virtual void destroy();

        // Legacy model: one lumped volume step at the fluid sub-step rate.
        // 1-D model: no-op (see processNetwork).
        void process(double dt);

        // 1-D model ---------------------------------------------------------
        inline bool isLegacy() const { return m_legacy; }
        // Register a cylinder's exhaust port; returns the primary index.
        // portArea is the head's exhaust runner cross-section, kFlow the legacy
        // primary flow constant (effective area = kFlow * sqrt(molar mass)).
        int attachPrimary(GasSystem *port, double length, double portArea, double kFlow);
        void buildNetwork();
        inline bool isNetworkBuilt() const { return m_network.isBuilt(); }
        // Per fluid sub-step, called by the cylinder in place of the lumped primary -> collector flow
        void exchangePort(int primaryIndex, double dt);
        // Once per simulation step
        void processNetwork(double dt);
        // Monopole source at the open end (Pa at 1 m), from the last processNetwork
        inline double getOutletSoundSource() const { return m_network.radiatedPressure(); }
        inline ExhaustNetwork *getNetwork() { return &m_network; }
        inline const ExhaustNetwork *getNetwork() const { return &m_network; }
        inline const Parameters &getParameters() const { return m_parameters; }
        // -------------------------------------------------------------------

        inline int getIndex() const { return m_index; }
        inline double getLength() const { return m_length; }
        inline double getFlow() const { return m_flow; }
        inline double getAudioVolume() const { return m_audioVolume; }
        inline double getPrimaryFlowRate() const { return m_primaryFlowRate; }
        inline double getCollectorCrossSectionArea() const { return m_collectorCrossSectionArea; }
        inline double getPrimaryTubeLength() const { return m_primaryTubeLength; }
        inline double getVelocityDecay() const { return m_velocityDecay; }
        inline ImpulseResponse *getImpulseResponse() const { return m_impulseResponse; }

        inline GasSystem *getSystem() { return &m_system; }

    protected:
        GasSystem m_atmosphere;
        GasSystem m_system;

        ImpulseResponse *m_impulseResponse;

        double m_length;
        double m_primaryTubeLength;
        double m_collectorCrossSectionArea;
        double m_primaryFlowRate;
        double m_outletFlowRate;
        double m_audioVolume;
        double m_velocityDecay;
        int m_index;

        double m_flow;

        bool m_legacy;
        Parameters m_parameters;
        ExhaustNetwork m_network;
};

#endif /* ATG_ENGINE_SIM_EXHAUST_SYSTEM_H */
