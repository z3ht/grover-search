#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <string>

#include "quantum_viz.hpp"

using Complex = std::complex<double>;
constexpr double PI = 3.1415926535897932384626433832795;

void apply_hadamard_on_qubit(std::vector<Complex>& psi, int q) {
    const int N = (int)psi.size();
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    const int stride = 1 << (q + 1);
    const int half   = 1 << q;
    for (int base = 0; base < N; base += stride) {
        for (int i = 0; i < half; ++i) {
            int i0 = base + i;
            int i1 = base + i + half;
            Complex a = psi[i0], b = psi[i1];
            psi[i0] = (a + b) * inv_sqrt2;
            psi[i1] = (a - b) * inv_sqrt2;
        }
    }
}

void apply_hadamard_all(std::vector<Complex>& psi, int n) {
    for (int q = 0; q < n; ++q) apply_hadamard_on_qubit(psi, q);
}

void apply_pauli_x_on_qubit(std::vector<Complex>& psi, int q) {
    const int N = (int)psi.size();
    const int mask = 1 << q;
    for (int i = 0; i < N; ++i) {
        int j = i ^ mask;
        if (i < j) std::swap(psi[i], psi[j]);
    }
}

template <class Pred>
void apply_oracle_phase_flip(std::vector<Complex>& psi, Pred f) {
    const int N = (int)psi.size();
    for (int i = 0; i < N; ++i) if (f(i)) psi[i] = -psi[i];
}

void apply_diffusion_mean(std::vector<Complex>& psi) {
    Complex mean = {0.0, 0.0};
    for (const auto& a : psi) mean += a;
    mean /= static_cast<double>(psi.size());
    for (auto& a : psi) a = 2.0 * mean - a;
}

int grover_iterations(int N, int M) {
    if (M <= 0) return 0;
    double k = std::floor((PI / 4.0) * std::sqrt(static_cast<double>(N) / static_cast<double>(M)));
    return static_cast<int>(std::max(0.0, k));
}

double grover_theta(int N, int M) {
    return std::asin(std::sqrt(static_cast<double>(M) / static_cast<double>(N)));
}
double grover_expected_marked_prob(int j, double theta) {
    double ang = (2 * j + 1) * theta;
    return std::sin(ang) * std::sin(ang);
}
double l2_from_uniform(const std::vector<Complex>& psi) {
    const double N = (double)psi.size();
    const Complex u = Complex(1.0 / std::sqrt(N), 0.0);
    double acc = 0.0;
    for (auto& a : psi) { Complex d = a - u; acc += std::norm(d); }
    return std::sqrt(acc);
}

struct Frame {
    std::vector<Complex> psi;
    std::vector<Complex> prev_for_reflection;
    std::string title_grid, title_arrows, title_bars, title_circuit;
    bool draw_reflection = false;
    bool show_circuit    = false;
    quantum_viz::OraclePhase oracle_phase = quantum_viz::OraclePhase::None;
    int  iter_index      = 0;
    int  highlight_index = -1; // which basis index to highlight this frame
};

static inline bool bit_is_zero_in_target(int target, int q) {
    return ((target >> q) & 1) == 0;
}

static void add_oracle_circuit_frames(std::vector<Frame>& frames,
                                      std::vector<Complex>& cur,
                                      int n, int target, int iter)
{
    std::vector<int> flip_bits;
    for (int q = 0; q < n; ++q)
        if (bit_is_zero_in_target(target, q)) flip_bits.push_back(q);

    // Target mapping under pre-X
    int transformed_target = target;
    for (int q : flip_bits) transformed_target ^= (1 << q);
    const int all_ones_state = (1 << n) - 1; // |11..1>

    // pre-X
    std::vector<Complex> before_prex = cur;
    for (int q : flip_bits) apply_pauli_x_on_qubit(cur, q);
    frames.push_back(Frame{
            cur, before_prex,
            "Oracle: pre-X (grid)",
            "Oracle: pre-X",
            "Oracle: pre-X (probs)",
            "Oracle step: pre-X",
            false, true, quantum_viz::OraclePhase::PreX, iter,
            transformed_target
    });

    // MCZ on |11..1>
    std::vector<Complex> before_mcz = cur;
    quantum_viz::apply_mcz_all_ones(cur, n);
    frames.push_back(Frame{
            cur, before_mcz,
            "Oracle: multi-control (grid)",
            "Oracle: multi-control",
            "Oracle: multi-control (probs)",
            "Oracle step: MCZ",
            false, true, quantum_viz::OraclePhase::MCX, iter,
            all_ones_state
    });

    // un-X
    std::vector<Complex> before_unx = cur;
    for (int q : flip_bits) apply_pauli_x_on_qubit(cur, q);
    frames.push_back(Frame{
            cur, before_unx,
            "Oracle: un-X (grid)",
            "Oracle: un-X",
            "Oracle: un-X (probs)",
            "Oracle step: un-X",
            false, true, quantum_viz::OraclePhase::UnX, iter,
            target
    });

    // summary frame after oracle
    frames.push_back(Frame{
            cur, {},
            "Oracle flip — iter " + std::to_string(iter) + " (grid)",
            "Oracle: phase flip — iter " + std::to_string(iter),
            "After oracle — iter " + std::to_string(iter) + " (probs)",
            "Oracle step: result",
            false, true, quantum_viz::OraclePhase::Result, iter,
            target
    });
}

int main() {
    const int n = 8;
    const int N = 1 << n;
    const int target = 5;

    auto is_marked = [target](int x) -> bool { return x == target; };
    int M = 0; for (int i = 0; i < N; ++i) if (is_marked(i)) ++M;
    if (M == 0) { std::cerr << "No marked items; exiting.\n"; return 0; }

    std::vector<Complex> psi(N, Complex{0.0, 0.0});
    psi[0] = Complex{1.0, 0.0};
    apply_hadamard_all(psi, n);

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

    const int k = grover_iterations(N, M);
    const double theta = grover_theta(N, M);

    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(6)
              << "N=" << N << ", M=" << M << ", k=" << k
              << ", theta=" << theta << " rad\n";

    std::vector<Frame> frames;
    frames.reserve(1 + (4 + 2) * k);

    frames.push_back(Frame{
            psi, {},
            "Init: uniform (grid)",
            "Init: uniform (Argand)",
            "Init: probabilities",
            "Oracle",
            false, false, quantum_viz::OraclePhase::None, 0,
            target
    });

    std::vector<Complex> cur = psi;
    for (int iter = 1; iter <= k; ++iter) {
        add_oracle_circuit_frames(frames, cur, n, target, iter);

        std::vector<Complex> after_oracle = cur;
        frames.push_back(Frame{
                after_oracle, after_oracle,
                "Diffusion (preview): reflect about mean",
                "Diffusion (preview): reflection target",
                "Diffusion (preview)",
                "Oracle",
                true, false, quantum_viz::OraclePhase::None, iter,
                target
        });

        apply_diffusion_mean(cur);
        frames.push_back(Frame{
                cur, after_oracle,
                "After diffusion — iter " + std::to_string(iter) + " (grid)",
                "After diffusion — iter " + std::to_string(iter),
                "After diffusion — iter " + std::to_string(iter) + " (probs)",
                "Oracle",
                false, false, quantum_viz::OraclePhase::None, iter,
                target
        });
    }

    const int F = (int)frames.size();
    int idx = 0;

    auto render = [&](int i){
        const Frame& f = frames[i];

        quantum_viz::GridArrowsConfig gcfg_copy = gcfg;
        gcfg_copy.marked_index = f.highlight_index;
        quantum_viz::show_grid_arrows(f.psi, n, gcfg_copy, f.title_grid.c_str());

        quantum_viz::ArrowsConfig acfg_copy = acfg;
        acfg_copy.marked_index = f.highlight_index;
        acfg_copy.draw_reflection = f.draw_reflection;
        quantum_viz::show_arrows(
                f.psi, n, acfg_copy, f.title_arrows.c_str(),
                f.draw_reflection ? &f.prev_for_reflection : nullptr
        );

        quantum_viz::ProbBarsConfig pcfg_copy = pcfg;
        pcfg_copy.marked_index = f.highlight_index;
        quantum_viz::show_prob_bars(f.psi, n, pcfg_copy, f.title_bars.c_str());

        if (f.show_circuit)
            quantum_viz::show_oracle_step(n, target, f.oracle_phase, ocfg,
                                          f.title_circuit.empty() ? "Oracle" : f.title_circuit.c_str());
        else
            quantum_viz::show_oracle_step(n, target, quantum_viz::OraclePhase::None, ocfg, " ");

        if (f.iter_index > 0) {
            const double p_marked = std::norm(f.psi[target]);
            const double p_any_theory = grover_expected_marked_prob(f.iter_index, theta);
            const double p_specific_theory = (M == 1) ? p_any_theory : p_any_theory / (double)M;
            const double l2 = l2_from_uniform(f.psi);
            std::cout << "Frame " << std::setw(2) << (i+1) << "/" << F
                      << "  Iter " << f.iter_index
                      << " | P[target]=" << p_marked
                      << " | theory≈ " << (M==1 ? p_specific_theory : p_any_theory)
                      << " | ||psi-|s>||2=" << l2 << "\n";
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
