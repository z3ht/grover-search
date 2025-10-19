#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include "quantum_viz.hpp"
#include "qcommon.hpp"
#include "frames.hpp"
#include "ialgorithm.hpp"
#include "shor_algorithm.hpp"
#include "grover_algorithm.hpp"


static quantum_viz::OraclePhase to_viz_phase(qframes::PhaseTag t) {
    using PT = qframes::PhaseTag;
    switch (t) {
        case PT::PreX:         return quantum_viz::OraclePhase::PreX;
        case PT::MCZ:          return quantum_viz::OraclePhase::MCX;
        case PT::UnX:          return quantum_viz::OraclePhase::UnX;
        case PT::OracleResult: return quantum_viz::OraclePhase::Result;
        case PT::None: default:return quantum_viz::OraclePhase::None;
    }
}

static int infer_n(const std::vector<qcommon::Complex>& psi) {
    int n = 0;
    while ((1 << n) < static_cast<int>(psi.size())) ++n;
    return n;
}

int main() {
    const int viz_qubits = 8;   // visualization grid size (2^n bins)
    const int target = 21;

    qcommon::QRegister reg(viz_qubits);
    std::unique_ptr<IAlgorithm> algo = std::make_unique<GroverAlgorithm>(target);

    auto frames = algo->run_and_make_frames(std::move(reg));

    if (frames.empty()) {
        std::cerr << "No frames produced.\n";
        return 0;
    }

    // Visualization configs
    quantum_viz::GridArrowsConfig gcfg;
    gcfg.r_bits = viz_qubits / 2;
    gcfg.cell = 28;
    gcfg.fill_prob = true;
    gcfg.hue_by_phase = false;
    gcfg.two_tone_sign = false;
    gcfg.gamma = 0.5;
    gcfg.min_frame_ms = 16;
    gcfg.arrow_len_mode = 0;

    quantum_viz::ArrowsConfig acfg;
    acfg.panel_size = 560;
    acfg.margin = 16;
    acfg.min_frame_ms = 16;
    acfg.draw_all = true;
    acfg.draw_mean = true;
    acfg.draw_reflection = false;
    acfg.draw_theta_arc = false;

    quantum_viz::ProbBarsConfig pcfg;
    pcfg.W = 560;
    pcfg.H = 300;
    pcfg.log_scale = false;
    pcfg.gamma = 0.55;
    pcfg.sort_descending = false;
    pcfg.min_frame_ms = 16;

    quantum_viz::OraclePanelConfig ocfg;
    ocfg.W = 720;
    ocfg.H = 240;
    ocfg.min_frame_ms = 16;
    ocfg.margin = 28;

    const int F = static_cast<int>(frames.size());
    int idx = 0;

    auto render = [&](int i) {
        const auto& f = frames[i];
        const int n_frame = infer_n(f.psi);

        auto gcfg_copy = gcfg;
        gcfg_copy.marked_index = f.highlight_index;
        quantum_viz::show_grid_arrows(f.psi, n_frame, gcfg_copy, f.title_grid.c_str());

        auto acfg_copy = acfg;
        acfg_copy.marked_index = f.highlight_index;
        acfg_copy.draw_reflection = f.draw_reflection;
        quantum_viz::show_arrows(
                f.psi, n_frame, acfg_copy, f.title_arrows.c_str(),
                f.draw_reflection ? &f.prev_for_reflection : nullptr
        );

        auto pcfg_copy = pcfg;
        pcfg_copy.marked_index = f.highlight_index;
        quantum_viz::show_prob_bars(f.psi, n_frame, pcfg_copy, f.title_bars.c_str());

        quantum_viz::BlochSpheresConfig bcfg;
        bcfg.sphere_size = 160;
        bcfg.cols = std::min(n_frame, 5);  // max 5 columns
        bcfg.min_frame_ms = 16;
        quantum_viz::show_bloch_spheres(f.psi, n_frame, bcfg, "Qubit Bloch Spheres");

        if (f.show_circuit)
            quantum_viz::show_oracle_step(n_frame, 0, to_viz_phase(f.phase), ocfg,
                                          f.title_circuit.empty() ? "Step" : f.title_circuit.c_str());
        else
            quantum_viz::show_oracle_step(n_frame, 0, quantum_viz::OraclePhase::None, ocfg, " ");

        std::cout << "Frame " << std::setw(2) << (i + 1) << "/" << F
                  << "  —  " << (f.title_grid.empty() ? "(frame)" : f.title_grid)
                  << "\n";
    };

    render(idx);
    for (;;) {
        auto step = quantum_viz::wait_for_step();
        if (step == quantum_viz::Step::Quit) break;
        if (step == quantum_viz::Step::Next)      { if (idx + 1 < F) { ++idx; render(idx); } }
        else if (step == quantum_viz::Step::Prev) { if (idx - 1 >= 0) { --idx; render(idx); } }
        else if (step == quantum_viz::Step::First){ idx = 0;         render(idx); }
        else if (step == quantum_viz::Step::Last) { idx = F - 1;     render(idx); }
    }

    quantum_viz::shutdown();
    return 0;
}
