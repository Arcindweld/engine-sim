#ifndef ATG_ENGINE_SIM_EXHAUST_SYSTEM_NODE_H
#define ATG_ENGINE_SIM_EXHAUST_SYSTEM_NODE_H

#include "object_reference_node.h"

#include "engine_context.h"
#include "impulse_response_node.h"

#include "engine_sim.h"

namespace es_script {

    class ExhaustSystemNode : public ObjectReferenceNode<ExhaustSystemNode> {
    public:
        ExhaustSystemNode() { /* void */ }
        virtual ~ExhaustSystemNode() { /* void */ }

        ExhaustSystem *generate(EngineContext *context) {
            ExhaustSystem *exhaust = context->getExhaust(this);
            ExhaustSystem::Parameters parameters = m_parameters;
            parameters.impulseResponse = m_impulseResponse->generate(context);

            exhaust->initialize(parameters);

            return exhaust;
        }

    protected:
        virtual void registerInputs() {
            addInput("length", &m_parameters.length);
            addInput("collector_cross_section_area", &m_parameters.collectorCrossSectionArea);
            addInput("outlet_flow_rate", &m_parameters.outletFlowRate);
            addInput("primary_tube_length", &m_parameters.primaryTubeLength);
            addInput("primary_flow_rate", &m_parameters.primaryFlowRate);
            addInput("audio_volume", &m_parameters.audioVolume);
            addInput("velocity_decay", &m_parameters.velocityDecay);
            addInput("impulse_response", &m_impulseResponse, InputTarget::Type::Object);

            // 1-D gas-dynamics model (optional in scripts; defaults in objects.mr)
            addInput("exhaust_model", &m_parameters.exhaustModel);
            addInput("legacy_exhaust", &m_parameters.legacy);
            addInput("primary_diameter", &m_parameters.primaryDiameter);
            addInput("collector_diameter", &m_parameters.collectorDiameter);
            addInput("collector_length", &m_parameters.collectorLength);
            addInput("megaphone_length", &m_parameters.megaphoneLength);
            addInput("megaphone_exit_diameter", &m_parameters.megaphoneExitDiameter);
            addInput("tailpipe_length", &m_parameters.tailpipeLength);
            addInput("tailpipe_diameter", &m_parameters.tailpipeDiameter);
            addInput("pipe_cell_size", &m_parameters.pipeCellSize);
            addInput("wall_temperature", &m_parameters.wallTemperature);

            ObjectReferenceNode<ExhaustSystemNode>::registerInputs();
        }

        virtual void _evaluate() {
            setOutput(this);

            // Read inputs
            readAllInputs();
        }

        ImpulseResponseNode *m_impulseResponse = nullptr;
        ExhaustSystem::Parameters m_parameters;
    };

} /* namespace es_script */

#endif /* ATG_ENGINE_SIM_EXHAUST_SYSTEM_NODE_H */
