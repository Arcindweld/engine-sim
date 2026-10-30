// Engine-level tests for the 1-D exhaust model, driven through the script
// compiler: the shipped engines load and run, legacy mode keeps the lumped
// model, primary length changes the breathing of the 1-D model only.
//
// Copyright (c) 2026 Arcindweld. MIT License (see LICENSE).
//
// The test executable needs the repository root (for es/ and assets/); it is
// passed by CMake as ENGINE_SIM_SOURCE_DIR and made the working directory.

#include <gtest/gtest.h>

#include "../scripting/include/compiler.h"
#include "../include/engine.h"
#include "../include/piston_engine_simulator.h"
#include "../include/exhaust_system.h"
#include "../include/units.h"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef ENGINE_SIM_SOURCE_DIR
#define ENGINE_SIM_SOURCE_DIR "."
#endif

namespace {

struct RunResult {
    bool compiled = false;
    bool nan = false;
    bool started = false;
    double meanTorque = 0.0;     // N m, dyno hold, last 0.5 s
    double meanRpm = 0.0;
    int cells = 0;
    double massImbalance = 0.0;  // relative to the throughput, 1-D systems only
    bool legacy = false;
};

// Compile a script, crank it for up to `crank` seconds at WOT, then hold it on
// the dyno at `holdRpm` for `hold` seconds.
RunResult runScript(const std::string &path, double crank, double holdRpm, double hold) {
    RunResult res;
    std::srand(12345);
    es_script::Compiler compiler;
    compiler.initialize();
    if (!compiler.compile(path.c_str())) {
        std::ifstream log("error_log.log");
        std::cerr << path << ": compile failed\n" << log.rdbuf() << "\n";
        compiler.destroy();
        return res;
    }
    const es_script::Compiler::Output output = compiler.execute();
    Engine *engine = output.engine;
    if (engine == nullptr || output.vehicle == nullptr || output.transmission == nullptr) {
        compiler.destroy();
        return res;
    }
    res.compiled = true;

    Simulator *sim = engine->createSimulator(output.vehicle, output.transmission);
    sim->setSimulationFrequency(static_cast<int>(engine->getSimulationFrequency()));
    engine->getIgnitionModule()->m_enabled = true;
    sim->m_starterMotor.m_enabled = true;
    for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
        ExhaustSystem *ex = engine->getExhaustSystem(e);
        res.legacy = res.legacy || ex->isLegacy();
        if (!ex->isLegacy()) res.cells += ex->getNetwork()->cellCount();
    }

    Synthesizer &syn = sim->synthesizer();
    const double startRpm = 2.5 * std::abs(units::toRpm(engine->getStarterSpeed()));
    double t = 0, speed = 0, torqueSum = 0, rpmSum = 0;
    int nAvg = 0;
    bool dynoOn = false;
    double tStart = -1;
    const double tEnd = crank + hold;
    while (t < tEnd) {
        if (!res.started && engine->getRpm() > startRpm) {
            res.started = true;
            tStart = t;
        }
        if (res.started || t > crank) sim->m_starterMotor.m_enabled = false;
        if (res.started && !dynoOn && holdRpm > 0) {
            sim->m_dyno.m_enabled = true;
            sim->m_dyno.m_hold = true;
            sim->m_dyno.m_rotationSpeed = units::rpm(holdRpm);
            dynoOn = true;
        }
        if (res.started && t > tStart + hold) break;
        speed = 0.5 + 0.5 * speed;
        engine->setSpeedControl(speed);
        sim->startFrame(1.0 / 60.0);
        while (sim->simulateStep()) { /* void */ }
        sim->endFrame();
        // drain the synthesizer input so startFrame() keeps its step count
        for (int c = 0; c < syn.m_inputChannelCount; ++c) {
            const size_t n = syn.m_inputChannels[c].data.size();
            if (n > 1000) syn.m_inputChannels[c].data.removeBeginning(n - 1000);
        }
        syn.m_latency = static_cast<int>(syn.m_inputChannels[0].data.size());
        t += sim->simulationSteps() * sim->getTimestep();

        const double rpm = engine->getRpm();
        const double torque = sim->getFilteredDynoTorque();
        if (std::isnan(rpm) || std::isnan(torque)) res.nan = true;
        for (int c = 0; c < engine->getCylinderCount(); ++c) {
            if (std::isnan(engine->getChamber(c)->m_system.pressure())) res.nan = true;
        }
        for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
            ExhaustSystem *ex = engine->getExhaustSystem(e);
            if (!ex->isLegacy() && std::isnan(ex->getNetwork()->totalMass())) res.nan = true;
        }
        if (res.started && t > tStart + hold - 0.5) {
            torqueSum += torque;
            rpmSum += rpm;
            ++nAvg;
        }
    }
    if (nAvg > 0) {
        res.meanTorque = torqueSum / nAvg;
        res.meanRpm = rpmSum / nAvg;
    }
    for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
        ExhaustSystem *ex = engine->getExhaustSystem(e);
        if (ex->isLegacy()) continue;
        const ExhaustNetwork *net = ex->getNetwork();
        const double imbalance = net->totalMass() - net->initialMass() - net->massIn() + net->massOut();
        if (net->massIn() != 0) res.massImbalance = std::max(res.massImbalance, std::fabs(imbalance / net->massIn()));
    }
    compiler.destroy();
    return res;
}

std::string writeScript(const std::string &name, const std::string &body) {
    const std::string path = name;
    std::ofstream(path) << body;
    return path;
}

std::string readFile(const std::string &path) {
    std::ifstream f(path);
    std::string s((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    return s;
}

// Insert `, <input>` after every impulse_response argument of a script.
std::string withExhaustInput(std::string s, const std::string &input) {
    size_t pos = 0;
    while ((pos = s.find("impulse_response:", pos)) != std::string::npos) {
        size_t end = pos + 17;
        while (end < s.size() && s[end] == ' ') ++end;
        int depth = 0;
        while (end < s.size()) {
            const char ch = s[end];
            if (ch == '(') ++depth;
            else if (ch == ')') { if (depth == 0) break; --depth; }
            else if ((ch == ',' || ch == '\n' || ch == '/') && depth == 0) break;
            ++end;
        }
        // trim trailing spaces of the expression
        size_t e = end;
        while (e > pos && s[e - 1] == ' ') --e;
        s.insert(e, ", " + input);
        pos = e + 2 + input.size();
    }
    return s;
}

std::string legacyVariant(const std::string &s) { return withExhaustInput(s, "legacy_exhaust: true"); }
std::string oneDVariant(const std::string &s) { return withExhaustInput(s, "exhaust_model: \"1d\""); }

struct ShippedEngine {
    const char *file;
    const char *node;   // nullptr: the file defines main()
};

const ShippedEngine shippedEngines[] = {
    { "atg-video-1/01_honda_trx520.mr", nullptr },
    { "atg-video-1/02_kohler_ch750.mr", nullptr },
    { "atg-video-1/03_harley_davidson_shovelhead.mr", nullptr },
    { "atg-video-1/04_hayabusa.mr", nullptr },
    { "atg-video-1/05_honda_vtec.mr", nullptr },
    { "atg-video-1/06_subaru_ej25.mr", nullptr },
    { "atg-video-1/07_audi_i5.mr", nullptr },
    { "atg-video-1/08_radial_5.mr", nullptr },
    { "atg-video-2/01_subaru_ej25_eh.mr", nullptr },
    { "atg-video-2/02_subaru_ej25_uh.mr", nullptr },
    { "atg-video-2/03_2jz.mr", nullptr },
    { "atg-video-2/04_60_degree_v6.mr", nullptr },
    { "atg-video-2/05_odd_fire_v6.mr", nullptr },
    { "atg-video-2/06_even_fire_v6.mr", nullptr },
    { "atg-video-2/07_gm_ls.mr", nullptr },
    { "atg-video-2/08_ferrari_f136_v8.mr", nullptr },
    { "atg-video-2/09_radial_9.mr", nullptr },
    { "atg-video-2/10_lfa_v10.mr", nullptr },
    { "atg-video-2/11_merlin_v12.mr", nullptr },
    { "atg-video-2/12_ferrari_412_t2.mr", nullptr },
    { "audi/i5.mr", "audi_i5_2_2L" },
    { "bmw/M52B28.mr", "M52B28" },
    { "chevrolet/chev_truck_454.mr", "chev_truck_454" },
    { "chevrolet/engine_03_for_e1.mr", "engine_03_for_e1" },
    { "kohler/kohler_ch750.mr", "kohler_ch750" },
};

std::string scriptFor(const ShippedEngine &e) {
    std::string s = "import \"engine_sim.mr\"\nimport \"assets/engines/" + std::string(e.file) + "\"\n";
    if (e.node == nullptr) s += "main()\n";
    else s += "set_engine(" + std::string(e.node) + "())\nset_vehicle(vehicle())\nset_transmission(transmission())\n";
    return s;
}

class EngineTests : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            std::filesystem::current_path(ENGINE_SIM_SOURCE_DIR);
        }
};

} // namespace

// (f) Every shipped engine, forced to the 1-D model (`exhaust_model: "1d"`
// inserted into a copy of its file), compiles and runs without NaN; the ones
// that start within the cranking time make positive torque on the dyno.  The
// Merlin V12 does not start in this headless start sequence with either
// exhaust model (it needs a longer crank than this test allows), so it is
// only checked for NaN.
TEST_F(EngineTests, ShippedEnginesRunWithOneDModel) {
    int started = 0, total = 0;
    for (const ShippedEngine &e : shippedEngines) {
        // the copy lives in the repository root: make the radials' local import absolute
        std::string text = oneDVariant(readFile(std::string("assets/engines/") + e.file));
        const std::string dir = std::string(e.file).substr(0, std::string(e.file).find('/'));
        for (size_t at; (at = text.find("import \"radial.mr\"")) != std::string::npos;) {
            text.replace(at, 18, "import \"assets/engines/" + dir + "/radial.mr\"");
        }
        std::ofstream("test_shipped_engine_tmp.mr") << text;
        std::string script = "import \"engine_sim.mr\"\nimport \"test_shipped_engine_tmp.mr\"\n";
        if (e.node == nullptr) script += "main()\n";
        else script += "set_engine(" + std::string(e.node) + "())\nset_vehicle(vehicle())\nset_transmission(transmission())\n";
        const std::string path = writeScript("test_shipped_tmp.mr", script);
        const RunResult r = runScript(path, 2.0, 2500.0, 1.5);
        ++total;
        std::cout << "[  values  ] " << e.file << ": cells " << r.cells
                  << (r.started ? ", started, torque " : ", not started in 2 s, torque ")
                  << r.meanTorque << " N m at " << r.meanRpm << " rpm, mass imbalance "
                  << r.massImbalance << "\n";
        EXPECT_TRUE(r.compiled) << e.file;
        EXPECT_FALSE(r.nan) << e.file;
        EXPECT_FALSE(r.legacy) << e.file << " was forced to the 1-D model";
        if (r.started) {
            ++started;
            EXPECT_GT(r.meanTorque, 0.0) << e.file;
            EXPECT_LT(r.massImbalance, 0.005) << e.file;
        }
    }
    std::remove("test_shipped_tmp.mr");
    std::remove("test_shipped_engine_tmp.mr");
    EXPECT_GE(started, total - 2);
}

namespace {
    // compile a script and report whether every exhaust system is legacy / 1-D
    void exhaustModes(const std::string &path, int *legacyCount, int *oneDCount) {
        *legacyCount = *oneDCount = 0;
        es_script::Compiler compiler;
        compiler.initialize();
        ASSERT_TRUE(compiler.compile(path.c_str())) << path;
        const es_script::Compiler::Output output = compiler.execute();
        ASSERT_NE(output.engine, nullptr);
        Engine *engine = output.engine;
        engine->createSimulator(output.vehicle, output.transmission);
        for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
            ExhaustSystem *ex = engine->getExhaustSystem(e);
            if (ex->isLegacy()) ++*legacyCount; else ++*oneDCount;
            EXPECT_EQ(ex->isNetworkBuilt(), !ex->isLegacy());
        }
        compiler.destroy();
    }
}

// exhaust_model "auto" (the default): a shipped engine, which gives none of
// the 1-D geometry inputs, runs the legacy model; the calculator export in
// test/data, which gives primary_diameter etc., runs the 1-D model; "1d",
// "legacy" and legacy_exhaust: true force the choice.
TEST_F(EngineTests, AutoModelSelection) {
    int legacy = 0, oneD = 0;

    writeScript("test_auto_tmp.mr", "import \"engine_sim.mr\"\nimport \"assets/engines/atg-video-1/06_subaru_ej25.mr\"\nmain()\n");
    exhaustModes("test_auto_tmp.mr", &legacy, &oneD);
    EXPECT_EQ(legacy, 2);
    EXPECT_EQ(oneD, 0);

    // (the calculator export is a complete script that calls main() itself)
    exhaustModes("test/data/exhaust_tuner_i4.mr", &legacy, &oneD);
    EXPECT_EQ(legacy, 0);
    EXPECT_EQ(oneD, 1);

    const std::string ej25 = readFile("assets/engines/atg-video-1/06_subaru_ej25.mr");
    const std::string i4 = readFile("test/data/exhaust_tuner_i4.mr");
    struct Case { std::string text; const char *tag; int legacy, oneD; };
    const Case cases[] = {
        { oneDVariant(ej25), "ej25 forced 1d", 0, 2 },
        { withExhaustInput(ej25, "exhaust_model: \"legacy\""), "ej25 explicit legacy", 2, 0 },
        { legacyVariant(ej25), "ej25 legacy_exhaust alias", 2, 0 },
        { withExhaustInput(i4, "exhaust_model: \"legacy\""), "i4 forced legacy", 1, 0 },
        { legacyVariant(i4), "i4 legacy_exhaust alias", 1, 0 },
        { withExhaustInput(ej25, "pipe_cell_size: 4 * units.cm"), "ej25 with pipe_cell_size", 0, 2 },
    };
    for (const Case &c : cases) {
        const bool complete = c.text.find("\nmain()") != std::string::npos;   // calls main() itself
        if (complete) {
            writeScript("test_auto_tmp.mr", c.text);
        }
        else {
            std::ofstream("test_auto_engine_tmp.mr") << c.text;
            writeScript("test_auto_tmp.mr", "import \"engine_sim.mr\"\nimport \"test_auto_engine_tmp.mr\"\nmain()\n");
        }
        exhaustModes("test_auto_tmp.mr", &legacy, &oneD);
        EXPECT_EQ(legacy, c.legacy) << c.tag;
        EXPECT_EQ(oneD, c.oneD) << c.tag;
    }
    std::remove("test_auto_tmp.mr");
    std::remove("test_auto_engine_tmp.mr");
}

// (g) legacy_exhaust: true keeps the lumped model: no network is built, the
// runner volume still holds the primary tube, the engine runs and makes
// torque.  (Bit-identity with the pre-change code was verified with two
// builds of the same harness, see docs/exhaust-1d.md.)
TEST_F(EngineTests, LegacyModeKeepsLumpedModel) {
    const std::string legacy = legacyVariant(readFile("assets/engines/atg-video-1/06_subaru_ej25.mr"));
    ASSERT_NE(legacy.find("legacy_exhaust: true"), std::string::npos);
    std::ofstream("test_legacy_engine_tmp.mr") << legacy;
    const std::string path = writeScript("test_legacy_tmp.mr",
        "import \"engine_sim.mr\"\nimport \"test_legacy_engine_tmp.mr\"\nmain()\n");

    // structural check
    {
        es_script::Compiler compiler;
        compiler.initialize();
        ASSERT_TRUE(compiler.compile(path.c_str()));
        const es_script::Compiler::Output output = compiler.execute();
        ASSERT_NE(output.engine, nullptr);
        Engine *engine = output.engine;
        Simulator *sim = engine->createSimulator(output.vehicle, output.transmission);
        (void)sim;
        for (int e = 0; e < engine->getExhaustSystemCount(); ++e) {
            EXPECT_TRUE(engine->getExhaustSystem(e)->isLegacy());
            EXPECT_FALSE(engine->getExhaustSystem(e)->isNetworkBuilt());
        }
        // runner volume = port volume + primary tube volume (combustion_chamber.cpp)
        CombustionChamber *ch = engine->getChamber(0);
        CylinderHead *head = ch->getCylinderHead();
        const int ci = ch->getPiston()->getCylinderIndex();
        ExhaustSystem *ex = head->getExhaustSystem(ci);
        const double expected = head->getExhaustRunnerVolume()
            + head->getExhaustRunnerCrossSectionArea() * (ex->getPrimaryTubeLength() + head->getHeaderPrimaryLength(ci));
        EXPECT_NEAR(ch->m_exhaustRunnerAndPrimary.volume(), expected, 1e-12);
        compiler.destroy();
    }

    const RunResult r = runScript(path, 2.0, 3000.0, 1.5);
    EXPECT_TRUE(r.compiled);
    EXPECT_FALSE(r.nan);
    EXPECT_TRUE(r.started);
    EXPECT_GT(r.meanTorque, 0.0);
    std::remove("test_legacy_tmp.mr");
    std::remove("test_legacy_engine_tmp.mr");
}

// Physics check: primary length changes the torque of the 1-D model (header
// tuning) but not that of the legacy model beyond its small volume effect.
// Uses the 2.0 L inline-4 with a 4-1 header (721 mm primaries, 44 mm) from
// test/data, held at 6000 rpm, with the primaries 10% shorter and 10% longer.
TEST_F(EngineTests, PrimaryLengthTunesOnlyTheOneDModel) {
    const std::string base = readFile("test/data/exhaust_tuner_i4.mr");
    const std::string key = "primary_tube_length: 720.969 * units.mm";
    ASSERT_NE(base.find(key), std::string::npos) << "test data changed; update the test";

    auto torqueFor = [&](const std::string &engineText, const char *tag) {
        const std::string path = std::string("test_tune_") + tag + ".mr";
        std::ofstream(path) << engineText;     // the file defines and calls main()
        const RunResult r = runScript(path, 2.0, 6000.0, 2.0);
        std::remove(path.c_str());
        EXPECT_TRUE(r.compiled && r.started && !r.nan) << tag;
        return r.meanTorque;
    };

    double tq1d[2], tqLegacy[2];
    const char *lengths[2] = { "primary_tube_length: 648.87 * units.mm", "primary_tube_length: 793.07 * units.mm" };
    for (int i = 0; i < 2; ++i) {
        std::string s = base;
        s.replace(s.find(key), key.size(), lengths[i]);
        tq1d[i] = torqueFor(s, i == 0 ? "short" : "long");
        tqLegacy[i] = torqueFor(legacyVariant(s), i == 0 ? "short_legacy" : "long_legacy");
    }
    const double change1d = std::fabs(tq1d[1] - tq1d[0]) / std::max(std::fabs(tq1d[0]), 1.0);
    const double changeLegacy = std::fabs(tqLegacy[1] - tqLegacy[0]) / std::max(std::fabs(tqLegacy[0]), 1.0);
    std::cout << "[  values  ] 6000 rpm, primaries -10% / +10%: 1-D " << tq1d[0] << " / " << tq1d[1]
              << " N m (" << 100 * change1d << " %), legacy " << tqLegacy[0] << " / " << tqLegacy[1]
              << " N m (" << 100 * changeLegacy << " %)\n";
    EXPECT_GT(change1d, 0.05);
    EXPECT_LT(changeLegacy, 0.03);
}
