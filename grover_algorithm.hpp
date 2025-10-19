#pragma once
#include "ialgorithm.hpp"
#include <iomanip>

struct GroverAlgorithm : public IAlgorithm {
    int target = 0;    // index in [0, 2^n - 1]
    int M = 1;        // number of marked items (default 1)

    explicit GroverAlgorithm(int target_idx, int marked_count = 1)
            : target(target_idx), M(marked_count) {}

    std::vector<qframes::Frame> run_and_make_frames(qcommon::QRegister reg) override {
        const int n = reg.n;
        const int N = static_cast<int>(reg.psi.size());
        if (target < 0 || target >= N) return {};

        // Init: |0...0> -> H^⊗n -> |s>
        reg.psi.assign(N, qcommon::Complex{0.0, 0.0});
        reg.psi[0] = qcommon::Complex{1.0, 0.0};
        qcommon::hadamard_all(reg.psi, n);

        // Iteration count & theory params
        const int k = qcommon::grover_iterations_for(N, M);

        std::vector<qframes::Frame> frames;
        frames.reserve(1 + (4 + 2) * k);

        // Initial frame
        frames.push_back(qframes::Frame{
                reg.psi, {},
                "Init: uniform (grid)",
                "Init: uniform (Argand)",
                "Init: probabilities",
                "Algorithm",
                false, false, qframes::PhaseTag::None, 0,
                target
        });

        // Work buffer for iterations
        auto cur = reg.psi;

        for (int iter = 1; iter <= k; ++iter) {
            // --- Oracle sandwich ---
            const auto map = qcommon::build_oracle_map(n, target);

            // pre-X
            auto before_prex = cur;
            for (int q : map.flip_bits) qcommon::pauli_x_on_qubit(cur, q);
            frames.push_back(qframes::Frame{
                    cur, before_prex,
                    "Oracle: pre-X (grid)",
                    "Oracle: pre-X",
                    "Oracle: pre-X (probs)",
                    "Oracle step: Pre-X",
                    false, true, qframes::PhaseTag::PreX, iter,
                    map.transformed_target
            });

            // MCZ
            auto before_mcz = cur;
            qcommon::mcz_all_ones(cur, n);
            frames.push_back(qframes::Frame{
                    cur, before_mcz,
                    "Oracle: multi-control (grid)",
                    "Oracle: multi-control",
                    "Oracle: multi-control (probs)",
                    "Oracle step: MCZ",
                    false, true, qframes::PhaseTag::MCZ, iter,
                    map.all_ones_index
            });

            // un-X
            auto before_unx = cur;
            for (int q : map.flip_bits) qcommon::pauli_x_on_qubit(cur, q);
            frames.push_back(qframes::Frame{
                    cur, before_unx,
                    "Oracle: un-X (grid)",
                    "Oracle: un-X",
                    "Oracle: un-X (probs)",
                    "Oracle step: Un-X",
                    false, true, qframes::PhaseTag::UnX, iter,
                    target
            });

            // summary/result
            frames.push_back(qframes::Frame{
                    cur, {},
                    "Oracle flip — iter " + std::to_string(iter) + " (grid)",
                    "Oracle: phase flip — iter " + std::to_string(iter),
                    "After oracle — iter " + std::to_string(iter) + " (probs)",
                    "Oracle step: result",
                    false, true, qframes::PhaseTag::OracleResult, iter,
                    target
            });

            // --- Diffusion ---
            auto after_oracle = cur;
            frames.push_back(qframes::Frame{
                    after_oracle, after_oracle,
                    "Diffusion (preview): reflect about mean",
                    "Diffusion (preview): reflection target",
                    "Diffusion (preview)",
                    "Algorithm",
                    true, false, qframes::PhaseTag::None, iter,
                    target
            });

            qcommon::diffusion_about_mean(cur);

            frames.push_back(qframes::Frame{
                    cur, after_oracle,
                    "After diffusion — iter " + std::to_string(iter) + " (grid)",
                    "After diffusion — iter " + std::to_string(iter),
                    "After diffusion — iter " + std::to_string(iter) + " (probs)",
                    "Algorithm",
                    false, false, qframes::PhaseTag::None, iter,
                    target
            });
        }

        return frames;
    }
};
