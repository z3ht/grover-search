#include <iostream>
#include <iomanip>
#include "quantum_viz.hpp"
#include "qcommon.hpp"
#include "frames.hpp"
#include "ialgorithm.hpp"
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

int main() {
    // ----- Configure problem -----
    const int n = 8;
    const int N = 1 << n;
    const int target = 5;       // try any in [0, N-1]

    // ----- Build initial register -----
    qcommon::QRegister reg(n);

    // ----- Choose algorithm (Grover for now; Shor later) -----
    std::unique_ptr<IAlgorithm> algo = std::make_unique<GroverAlgorithm>(target /*, M=1*/);

    // ----- Run algorithm and get frames (no viz performed here) -----
    std::vector<qframes::Frame> frames = algo->run_and_make_frames(std::move(reg));

    if (frames.empty()) {
        std::cerr << "No frames produced.\n";
        return 0;
    }

    // ----- Viz configs (owned only by this file) -----
    quantum_viz::GridArrowsConfig gcfg;
    gcfg.r_bits = n / 2; gcfg.cell = 28; gcfg.marked_index = target;
    gcfg.fill_prob = true; gcfg.hue_by_phase = false; gcfg.two_tone_sign = false;
    gcfg.gamma = 0.5; gcfg.min_frame_ms = 16; gcfg.arrow_len_mode = 0;

    quantum_viz::ArrowsConfig acfg;
    acfg.panel_size = 560; acfg.margin = 16; acfg.min_frame_ms = 16;
    acfg.marked_index = target; acfg.draw_all = true; acfg.draw_mean = true;
    acfg.draw_reflection = false; acfg.draw_theta_arc = false;

    quantum_viz::ProbBarsConfig pcfg;
    pcfg.W = 560; pcfg.H = 300; pcfg.log_scale = false; pcfg.gamma = 0.55;
    pcfg.sort_descending = false; pcfg.marked_index = target; pcfg.min_frame_ms = 16;

    quantum_viz::OraclePanelConfig ocfg;
    ocfg.W = 720; ocfg.H = 240; ocfg.min_frame_ms = 16; ocfg.margin = 28;

    // ----- Render loop -----
    const int F = static_cast<int>(frames.size());
    int idx = 0;

    auto render = [&](int i){
        const auto& f = frames[i];

        // grid
        auto gcfg_copy = gcfg;
        gcfg_copy.marked_index = f.highlight_index;
        quantum_viz::show_grid_arrows(f.psi, n, gcfg_copy, f.title_grid.c_str());

        // arrows
        auto acfg_copy = acfg;
        acfg_copy.marked_index = f.highlight_index;
        acfg_copy.draw_reflection = f.draw_reflection;
        quantum_viz::show_arrows(
                f.psi, n, acfg_copy, f.title_arrows.c_str(),
                f.draw_reflection ? &f.prev_for_reflection : nullptr
        );

        // bars
        auto pcfg_copy = pcfg;
        pcfg_copy.marked_index = f.highlight_index;
        quantum_viz::show_prob_bars(f.psi, n, pcfg_copy, f.title_bars.c_str());

        // circuit panel (optional)
        if (f.show_circuit)
            quantum_viz::show_oracle_step(n, target, to_viz_phase(f.phase), ocfg,
                                          f.title_circuit.empty() ? "Oracle" : f.title_circuit.c_str());
        else
            quantum_viz::show_oracle_step(n, target, quantum_viz::OraclePhase::None, ocfg, " ");

        // console readout (uses original target)
        if (f.iter_index > 0) {
            const double p_marked = std::norm(f.psi[target]);
            const double theta = qcommon::grover_theta_for(N, /*M=*/1);
            const double p_any_theory = qcommon::expected_marked_prob(f.iter_index, theta);
            std::cout << "Frame " << std::setw(2) << (i+1) << "/" << F
                      << "  Iter " << f.iter_index
                      << " | P[target]=" << p_marked
                      << " | theory≈ " << p_any_theory
                      << "\n";
        } else {
            std::cout << "Frame " << std::setw(2) << (i+1) << "/" << F << "  (init)\n";
        }
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
