// Headless Engine Simulator runner (Arcindweld fork, 1-D exhaust): no window, no audio device.
//
// Must run with the current directory containing an `es/` folder.
//
// usage: engine-sim-headless <file.mr> [--json out.json] [--sim seconds] [--throttle 0..1]
//        [--hold-rpm N] [--wav out.wav] [--no-audio] [--legacy] [--settle s]
//        [--cell-size m] [--sim-freq hz] [--trace out.csv] [--gain g]

#include "scripting/include/compiler.h"
#include "include/engine.h"
#include "include/piston_engine_simulator.h"
#include "include/exhaust_system.h"
#include "include/intake.h"
#include "include/cylinder_head.h"
#include "include/cylinder_bank.h"
#include "include/piston.h"
#include "include/ignition_module.h"
#include "include/units.h"
#include "include/constants.h"

#include <chrono>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::string jnum(double v) {
    if (!std::isfinite(v)) return "null";
    char b[64];
    std::snprintf(b, sizeof(b), "%.9g", v);
    return b;
}

bool readWav16(const std::string &path, std::vector<int16_t> *out) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    char riff[12];
    f.read(riff, 12);
    if (std::strncmp(riff, "RIFF", 4) != 0 || std::strncmp(riff + 8, "WAVE", 4) != 0) return false;
    uint16_t channels = 1, bits = 16;
    while (f) {
        char id[4]; uint32_t size = 0;
        f.read(id, 4); f.read(reinterpret_cast<char *>(&size), 4);
        if (!f) break;
        if (std::strncmp(id, "fmt ", 4) == 0) {
            std::vector<char> b(size);
            f.read(b.data(), size);
            std::memcpy(&channels, b.data() + 2, 2);
            std::memcpy(&bits, b.data() + 14, 2);
        }
        else if (std::strncmp(id, "data", 4) == 0) {
            if (bits != 16) return false;
            std::vector<int16_t> raw(size / 2);
            f.read(reinterpret_cast<char *>(raw.data()), size);
            out->clear();
            for (size_t i = 0; i < raw.size(); i += channels) out->push_back(raw[i]);
            return true;
        }
        else {
            f.seekg(size + (size & 1), std::ios::cur);
        }
    }
    return false;
}

bool writeWav16(const std::string &path, const std::vector<int16_t> &s, uint32_t rate) {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    const uint32_t dataBytes = static_cast<uint32_t>(s.size() * 2);
    const uint32_t riffSize = 36 + dataBytes;
    const uint16_t fmt = 1, ch = 1, bits = 16, align = 2;
    const uint32_t fmtSize = 16, byteRate = rate * 2;
    f.write("RIFF", 4); f.write(reinterpret_cast<const char *>(&riffSize), 4);
    f.write("WAVE", 4); f.write("fmt ", 4);
    f.write(reinterpret_cast<const char *>(&fmtSize), 4);
    f.write(reinterpret_cast<const char *>(&fmt), 2);
    f.write(reinterpret_cast<const char *>(&ch), 2);
    f.write(reinterpret_cast<const char *>(&rate), 4);
    f.write(reinterpret_cast<const char *>(&byteRate), 4);
    f.write(reinterpret_cast<const char *>(&align), 2);
    f.write(reinterpret_cast<const char *>(&bits), 2);
    f.write("data", 4); f.write(reinterpret_cast<const char *>(&dataBytes), 4);
    f.write(reinterpret_cast<const char *>(s.data()), dataBytes);
    return true;
}

// Access to protected members for the legacy override
struct ExhaustAccess : public ExhaustSystem {
    static void forceLegacy(ExhaustSystem *e, bool legacy) {
        ExhaustAccess *a = static_cast<ExhaustAccess *>(e);
        a->m_legacy = legacy;
        a->m_parameters.legacy = legacy;
    }
    static void setCellSize(ExhaustSystem *e, double cell) {
        static_cast<ExhaustAccess *>(e)->m_parameters.pipeCellSize = cell;
    }
};

CylinderHead *m_headOf(Engine *engine, int c) {
    Piston *piston = engine->getPiston(c);
    return engine->getHead(piston->getCylinderBank()->getIndex());
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "usage: engine-sim-headless <file.mr> [options]\n";
        return 2;
    }
    std::string script = argv[1], jsonPath, wavPath, tracePath;
    double simSeconds = 0, throttle = 1.0, holdRpm = 0, settle = 2.0, cellSize = 0, simFreq = 0, starterTime = 2.0, gain = 1.0;
    bool audio = true, forceLegacy = false, debugNan = false;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(i + 1 < argc ? argv[++i] : ""); };
        if (a == "--json") jsonPath = next();
        else if (a == "--wav") wavPath = next();
        else if (a == "--trace") tracePath = next();
        else if (a == "--sim") simSeconds = std::atof(next().c_str());
        else if (a == "--throttle") throttle = std::atof(next().c_str());
        else if (a == "--hold-rpm") holdRpm = std::atof(next().c_str());
        else if (a == "--settle") settle = std::atof(next().c_str());
        else if (a == "--cell-size") cellSize = std::atof(next().c_str());
        else if (a == "--sim-freq") simFreq = std::atof(next().c_str());
        else if (a == "--starter-time") starterTime = std::atof(next().c_str());
        else if (a == "--gain") gain = std::atof(next().c_str());   // synthesizer output volume (1 = engine-sim default; 0.4 avoids clipping)
        else if (a == "--no-audio") audio = false;
        else if (a == "--legacy") forceLegacy = true;
        else if (a == "--debug-nan") debugNan = true;
    }

    std::srand(12345);

    es_script::Compiler compiler;
    compiler.initialize();
    const bool compiled = compiler.compile(script.c_str());
    if (!compiled) {
        std::cout << "COMPILE: FAILED\n";
        std::ifstream log("error_log.log");
        std::cout << log.rdbuf() << "\n";
        compiler.destroy();
        return 1;
    }
    std::cout << "COMPILE: OK\n";
    const es_script::Compiler::Output output = compiler.execute();
    Engine *engine = output.engine;
    Vehicle *vehicle = output.vehicle;
    Transmission *transmission = output.transmission;
    std::cout << "EXECUTE: engine=" << (engine ? "yes" : "NO")
              << " vehicle=" << (vehicle ? "yes" : "NO")
              << " transmission=" << (transmission ? "yes" : "NO") << "\n";
    if (engine == nullptr || vehicle == nullptr || transmission == nullptr) {
        compiler.destroy();
        return 1;
    }

    // NOTE: --legacy here only switches the per-step processing and the sound
    // path; the chambers were already initialised by the script (runner
    // volumes, primary registration).  Use the script's legacy_exhaust input
    // for a faithful legacy run.  --cell-size takes effect because the
    // networks are built in createSimulator().
    for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
        if (cellSize > 0) ExhaustAccess::setCellSize(engine->getExhaustSystem(e), cellSize);
    }

    Simulator *sim = engine->createSimulator(vehicle, transmission);
    engine->calculateDisplacement();

    std::ostringstream js;
    js << "{\n  \"name\": \"" << engine->getName() << "\",\n";
    js << "  \"cylinders\": " << engine->getCylinderCount()
       << ", \"exhaust_systems\": " << engine->getExhaustSystemCount() << ",\n";
    js << "  \"displacement_L\": " << jnum(engine->getDisplacement() / units::L) << ",\n";
    js << "  \"sim_freq\": " << jnum(engine->getSimulationFrequency()) << ",\n";
    js << "  \"exhaust\": [";
    int totalCells = 0;
    for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
        ExhaustSystem *ex = engine->getExhaustSystem(e);
        js << (e ? ", " : "") << "{\"legacy\": " << (ex->isLegacy() ? "true" : "false")
           << ", \"length_mm\": " << jnum(ex->getLength() / units::mm)
           << ", \"collector_area_cm2\": " << jnum(ex->getCollectorCrossSectionArea() / (units::cm * units::cm))
           << ", \"primary_tube_length_mm\": " << jnum(ex->getPrimaryTubeLength() / units::mm)
           << ", \"audio_volume\": " << jnum(ex->getAudioVolume());
        if (!ex->isLegacy()) {
            const ExhaustNetwork *net = ex->getNetwork();
            totalCells += net->cellCount();
            js << ", \"primaries\": " << net->primaryCount()
               << ", \"collector_pipes\": " << net->collectorPipeCount()
               << ", \"cells\": " << net->cellCount()
               << ", \"cell_size_mm\": " << jnum(net->primaryPipe(0).cellSize() / units::mm)
               << ", \"primary_lengths_mm\": [";
            for (int k = 0; k < net->primaryCount(); ++k) {
                js << (k ? ", " : "") << jnum(net->primaryPipe(k).length() / units::mm);
            }
            js << "], \"primary_diameter_mm\": " << jnum(net->primaryPipe(0).faceArea(0) > 0 ? std::sqrt(4 * net->primaryPipe(0).faceArea(0) / constants::pi) / units::mm : 0)
               << ", \"collector_length_mm\": " << jnum(net->collectorPipe(0).length() / units::mm);
        }
        js << "}";
    }
    js << "],\n  \"total_cells\": " << totalCells;

    if (simSeconds > 0) {
        sim->setSimulationFrequency(static_cast<int>(simFreq > 0 ? simFreq : engine->getSimulationFrequency()));
        Synthesizer::AudioParameters ap = sim->synthesizer().getAudioParameters();
        ap.inputSampleNoise = static_cast<float>(engine->getInitialJitter());
        ap.airNoise = static_cast<float>(engine->getInitialNoise());
        ap.dF_F_mix = static_cast<float>(engine->getInitialHighFrequencyGain());
        ap.volume = static_cast<float>(gain);
        sim->synthesizer().setAudioParameters(ap);
        for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
            ImpulseResponse *ir = engine->getExhaustSystem(e)->getImpulseResponse();
            std::vector<int16_t> w;
            if (!readWav16(ir->getFilename(), &w)) {
                std::cerr << "IR not readable: " << ir->getFilename() << "\n";
                w.assign(1, INT16_MAX);
            }
            sim->synthesizer().initializeImpulseResponse(w.data(), (unsigned)w.size(), (float)ir->getVolume(), e);
        }

        IgnitionModule *im = engine->getIgnitionModule();
        im->m_enabled = true;
        sim->m_starterMotor.m_enabled = true;
        // The dyno picks its direction from the crank's current rotation, so
        // it is enabled only once the engine has started (see the loop).
        bool dynoArmed = holdRpm > 0;

        std::ofstream trace;
        if (!tracePath.empty()) {
            trace.open(tracePath);
            trace << "t,rpm,torque";
            for (int e = 0; e < engine->getExhaustSystemCount(); ++e) trace << ",p_out" << e << ",mdot_out" << e << ",src" << e;
            trace << "\n";
        }

        // Mass accounting for the 1-D networks: pipe mass + mass out - mass in - initial = 0
        std::vector<int16_t> audioSamples;
        const double frame = 1.0 / 60.0;
        double t = 0, speed = 0, sumRpm = 0, sumTorque = 0, minTorque = 1e30, maxTorque = -1e30;
        long long nAvg = 0;
        bool nanSeen = false;
        double simWall = 0, audioWall = 0, simCpu = 0;
        long long steps = 0;
        int subStepsMax = 0; double subStepsSum = 0; long long subStepsN = 0;
        Synthesizer &syn = sim->synthesizer();
        while (t < simSeconds) {
            if (t > starterTime || engine->getRpm() > 2.5 * std::abs(units::toRpm(engine->getStarterSpeed()))) {
                sim->m_starterMotor.m_enabled = false;
                if (dynoArmed && engine->getRpm() > 2.5 * std::abs(units::toRpm(engine->getStarterSpeed()))) {
                    sim->m_dyno.m_enabled = true;
                    sim->m_dyno.m_hold = true;
                    sim->m_dyno.m_rotationSpeed = units::rpm(holdRpm);
                    dynoArmed = false;
                }
            }
            speed = 0.5 * throttle + 0.5 * speed;
            engine->setSpeedControl(speed);

            auto s0 = std::chrono::steady_clock::now();
            const std::clock_t c0 = std::clock();
            sim->startFrame(frame);
            while (sim->simulateStep()) {
                ++steps;
                if (debugNan && !nanSeen) {
                    for (int c = 0; c < engine->getCylinderCount(); ++c) {
                        CombustionChamber *ch = engine->getChamber(c);
                        GasSystem &port = ch->m_exhaustRunnerAndPrimary;
                        bool pipeBad = false;
                        {
                            Piston *piston = engine->getPiston(c);
                            CylinderHead *head = m_headOf(engine, c);
                            ExhaustSystem *ex = head->getExhaustSystem(piston->getCylinderIndex());
                            if (!ex->isLegacy()) {
                                const ExhaustNetwork *net = ex->getNetwork();
                                for (int k = 0; k < net->primaryCount(); ++k) {
                                    const Pipe1D &P = net->primaryPipe(k);
                                    if (P.pressure(0) < 5000 || P.temperature(0) < 150 || P.temperature(0) > 4000) pipeBad = true;
                                }
                            }
                        }
                        if (pipeBad || std::isnan(ch->m_system.pressure()) || std::isnan(port.pressure()) || port.n() <= 0 || port.pressure() < 5000 || port.temperature() > 4000 || port.temperature() < 150) {
                            Piston *piston = engine->getPiston(c);
                            CylinderHead *head = m_headOf(engine, c);
                            ExhaustSystem *ex = head->getExhaustSystem(piston->getCylinderIndex());
                            std::cerr << "step " << steps << " t=" << t << " cyl " << c << ": cyl p=" << ch->m_system.pressure() << " T=" << ch->m_system.temperature()
                                      << " port p=" << port.pressure() << " n=" << port.n() << " T=" << port.temperature() << " V=" << port.volume() << "\n";
                            if (!ex->isLegacy()) {
                                const ExhaustNetwork *net = ex->getNetwork();
                                for (int k = 0; k < net->primaryCount(); ++k) {
                                    const Pipe1D &P = net->primaryPipe(k);
                                    std::cerr << "   primary " << k << " cell0 p=" << P.pressure(0) << " u=" << P.velocity(0) << " T=" << P.temperature(0) << " rho=" << P.density(0) << "\n";
                                }
                            }
                            nanSeen = true;
                            break;
                        }
                    }
                }
                if (trace.is_open() && steps % 10 == 0) {
                    trace << jnum(t + (steps % 100000) * 0) << "," << jnum(engine->getRpm()) << "," << jnum(sim->m_dyno.getTorque());
                    for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
                        ExhaustSystem *ex = engine->getExhaustSystem(e);
                        if (ex->isLegacy()) trace << "," << jnum(ex->getSystem()->pressure()) << ",0,0";
                        else trace << "," << jnum(ex->getNetwork()->outletPressure()) << "," << jnum(ex->getNetwork()->outletMassFlow()) << "," << jnum(ex->getOutletSoundSource());
                    }
                    trace << "\n";
                }
            }
            sim->endFrame();
            auto s1 = std::chrono::steady_clock::now();
            simWall += std::chrono::duration<double>(s1 - s0).count();
            simCpu += double(std::clock() - c0) / CLOCKS_PER_SEC;
            for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
                ExhaustSystem *ex = engine->getExhaustSystem(e);
                if (!ex->isLegacy()) {
                    const int ss = ex->getNetwork()->lastSubSteps();
                    subStepsMax = std::max(subStepsMax, ss); subStepsSum += ss; ++subStepsN;
                }
            }

            if (audio) {
                auto a0 = std::chrono::steady_clock::now();
                if (syn.m_inputChannels[0].data.size() > 0 && syn.m_audioBuffer.size() < 2000 && !syn.m_processed) {
                    syn.renderAudio();
                }
                int16_t buf[4096];
                const int n = sim->readAudioOutput(std::min<int>(4096, (int)syn.m_audioBuffer.size()), buf);
                audioSamples.insert(audioSamples.end(), buf, buf + n);
                auto a1 = std::chrono::steady_clock::now();
                audioWall += std::chrono::duration<double>(a1 - a0).count();
            }
            else {
                // keep the synthesizer input latency low so startFrame() does
                // not throttle the step count (it targets 0.1 s of latency)
                for (int c = 0; c < syn.m_inputChannelCount; ++c) {
                    const size_t n = syn.m_inputChannels[c].data.size();
                    if (n > 1000) syn.m_inputChannels[c].data.removeBeginning(n - 1000);
                }
                syn.m_latency = (int)syn.m_inputChannels[0].data.size();
            }

            t += sim->simulationSteps() * sim->getTimestep();
            const double rpm = engine->getRpm();
            const double torque = sim->getFilteredDynoTorque();
            if (std::isnan(rpm) || std::isnan(torque)) nanSeen = true;
            for (int c = 0; c < engine->getCylinderCount(); ++c) {
                CombustionChamber *ch = engine->getChamber(c);
                const bool bad = std::isnan(ch->m_system.pressure()) || std::isnan(ch->m_exhaustRunnerAndPrimary.pressure())
                    || ch->m_exhaustRunnerAndPrimary.n() <= 0 || ch->m_system.n() <= 0;
                if (bad && !nanSeen) {
                    std::cerr << "NaN/empty at t=" << t << " cyl " << c << " cyl p=" << ch->m_system.pressure()
                              << " n=" << ch->m_system.n() << " port p=" << ch->m_exhaustRunnerAndPrimary.pressure()
                              << " port n=" << ch->m_exhaustRunnerAndPrimary.n() << " T=" << ch->m_exhaustRunnerAndPrimary.temperature() << "\n";
                }
                if (bad) nanSeen = true;
            }
            for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
                ExhaustSystem *ex = engine->getExhaustSystem(e);
                if (!ex->isLegacy() && std::isnan(ex->getNetwork()->totalMass()) && !nanSeen) {
                    std::cerr << "network NaN at t=" << t << " exhaust " << e << "\n";
                    nanSeen = true;
                }
            }
            if (t > settle) {
                sumRpm += rpm; sumTorque += torque; ++nAvg;
                minTorque = std::min(minTorque, torque); maxTorque = std::max(maxTorque, torque);
            }
        }
        double peak = 0, rms = 0;
        for (int16_t s : audioSamples) { peak = std::max(peak, std::abs((double)s)); rms += (double)s * s; }
        rms = audioSamples.empty() ? 0 : std::sqrt(rms / audioSamples.size());

        js << ",\n  \"sim\": {\"seconds\": " << jnum(t) << ", \"throttle\": " << jnum(throttle)
           << ", \"hold_rpm\": " << jnum(holdRpm) << ", \"settle\": " << jnum(settle)
           << ", \"steps\": " << steps
           << ", \"mean_rpm\": " << jnum(nAvg ? sumRpm / nAvg : 0)
           << ", \"mean_dyno_torque_Nm\": " << jnum(nAvg ? sumTorque / nAvg : 0)
           << ", \"min_dyno_torque_Nm\": " << jnum(nAvg ? minTorque : 0)
           << ", \"max_dyno_torque_Nm\": " << jnum(nAvg ? maxTorque : 0)
           << ", \"nan\": " << (nanSeen ? "true" : "false")
           << ", \"peak_T_cyl0_K\": " << jnum(engine->getChamber(0)->m_peakTemperature)
           << ", \"fuel_consumed_g\": " << jnum(engine->getTotalFuelMassConsumed() / units::g)
           << ", \"manifold_kPa\": " << jnum(engine->getManifoldPressure() / units::kPa)
           << ", \"afr\": " << jnum(engine->getIntakeAfr())
           << ", \"sim_wall_s\": " << jnum(simWall) << ", \"audio_wall_s\": " << jnum(audioWall)
           << ", \"sim_cost_per_sim_second\": " << jnum(t > 0 ? simWall / t : 0)
           << ", \"sim_cpu_per_sim_second\": " << jnum(t > 0 ? simCpu / t : 0)
           << ", \"substeps_mean\": " << jnum(subStepsN ? subStepsSum / subStepsN : 0)
           << ", \"substeps_max\": " << subStepsMax
           << ", \"audio_samples\": " << audioSamples.size() << ", \"audio_peak\": " << jnum(peak)
           << ", \"audio_rms\": " << jnum(rms) << ",\n   \"mass\": [";
        for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
            ExhaustSystem *ex = engine->getExhaustSystem(e);
            js << (e ? ", " : "");
            if (ex->isLegacy()) { js << "null"; continue; }
            const ExhaustNetwork *net = ex->getNetwork();
            const double inPipes = net->totalMass(), in = net->massIn(), out = net->massOut(), init = net->initialMass();
            js << "{\"initial_kg\": " << jnum(init) << ", \"in_kg\": " << jnum(in) << ", \"out_kg\": " << jnum(out)
               << ", \"in_pipes_kg\": " << jnum(inPipes)
               << ", \"error_kg\": " << jnum(inPipes - init - in + out)
               << ", \"error_rel_to_throughput\": " << jnum((in != 0) ? (inPipes - init - in + out) / in : 0) << "}";
        }
        js << "]}";
        if (!wavPath.empty()) writeWav16(wavPath, audioSamples, 44100);
    }
    js << "\n}\n";

    if (!jsonPath.empty()) { std::ofstream(jsonPath) << js.str(); }
    else std::cout << js.str();

    compiler.destroy();
    return 0;
}
