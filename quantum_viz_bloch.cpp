// quantum_viz_bloch.cpp
#include "quantum_viz.hpp"
#include "bloch_viz.hpp"
#include <SDL.h>
#include <cmath>
#include <algorithm>

namespace quantum_viz {

    bool show_bloch_spheres(const std::vector<Complex>& psi, int n,
                            const BlochSpheresConfig& cfg, const char* title)
    {
        // Use the computation from bloch_viz.hpp
        auto qubits = qviz_local::compute_qubit_bloch_vectors(psi, n);

        const int rows = (n + cfg.cols - 1) / cfg.cols;
        const int W = cfg.cols * cfg.sphere_size + (cfg.cols + 1) * cfg.margin;
        const int H = rows * (cfg.sphere_size + 60) + cfg.margin * 2;

        if (!sdlm::ensure_window(sdlm::M().bloch, W, H, title)) return false;
        SDL_Surface* surf = sdlm::M().bloch.surf;
        SDL_FillRect(surf, nullptr, rgb(surf->format, 12, 12, 14));

        for (int q = 0; q < n; ++q) {
            const auto& v = qubits[q];

            int col = q % cfg.cols;
            int row = q / cfg.cols;
            int base_x = cfg.margin + col * (cfg.sphere_size + cfg.margin);
            int base_y = cfg.margin + row * (cfg.sphere_size + 60);

            int cx = base_x + cfg.sphere_size / 2;
            int cy = base_y + cfg.sphere_size / 2;
            int R = (cfg.sphere_size / 2) - 10;

            // Draw background circle (equator)
            circle(surf, cx, cy, R, rgb(surf->format, 60, 60, 70));

            if (cfg.draw_axes) {
                // X axis (horizontal)
                Uint32 axis_col = rgb(surf->format, 80, 80, 90);
                line(surf, cx - R, cy, cx + R, cy, axis_col);
                // Y axis (vertical, into/out of screen - show as dots)
                for (int i = -R; i <= R; i += 8) {
                    pset(surf, cx, cy + i, axis_col);
                }
                // Z axis (vertical)
                line(surf, cx, cy - R, cx, cy + R, axis_col);

                // Axis labels (simple dots for |0> and |1>)
                Uint32 label_col = rgb(surf->format, 150, 150, 160);
                // |0> at top
                boxfill(surf, cx - 3, cy - R - 8, 6, 6, label_col);
                // |1> at bottom
                boxfill(surf, cx - 3, cy + R + 2, 6, 6, label_col);
            }

            // Draw the Bloch vector as an arrow
            if (cfg.draw_state_vector && v.length > 1e-6) {
                // Convert spherical to 2D projection (XZ plane primarily)
                // theta: 0 = north pole |0>, pi = south pole |1>
                // phi: rotation around Z axis

                double x_3d = v.length * std::sin(v.theta) * std::cos(v.phi);
                double y_3d = v.length * std::sin(v.theta) * std::sin(v.phi);
                double z_3d = v.length * std::cos(v.theta);

                // Project onto 2D (simple orthographic: X horizontal, Z vertical, Y depth)
                int tip_x = cx + (int)(x_3d * R);
                int tip_y = cy - (int)(z_3d * R); // negative because screen Y is down

                // Color based on purity (length of Bloch vector)
                Uint32 vec_col = rgb(surf->format, 50, 200, 255);
                Uint32 head_col = rgb(surf->format, 80, 220, 255);

                arrow(surf, cx, cy, tip_x, tip_y, vec_col, head_col, 8.0, 22.0);

                // Draw a faint circle to show Y component (phase)
                if (std::abs(y_3d) > 0.01) {
                    int y_offset = (int)(y_3d * R * 0.3); // scale down for visibility
                    Uint32 phase_col = rgb(surf->format, 100, 100, 120);
                    circle(surf, cx + y_offset, cy, (int)(R * 0.95), phase_col);
                }
            }

            // Draw probability bar and info below sphere
            if (cfg.draw_probabilities) {
                int info_y = base_y + cfg.sphere_size + 10;

                // Qubit label
                Uint32 text_col = rgb(surf->format, 200, 200, 220);
                int label_x = base_x + 5;
                boxfill(surf, label_x, info_y, 30, 12, rgb(surf->format, 40, 40, 50));

                // P(|0>) and P(|1>) bars
                int bar_y = info_y + 18;
                int bar_w = cfg.sphere_size - 10;
                int bar_h = 8;

                // Background
                boxfill(surf, base_x + 5, bar_y, bar_w, bar_h, rgb(surf->format, 40, 40, 50));

                // P(|0>) in green-ish, P(|1>) in red-ish
                int p0_w = (int)(v.p0 * bar_w);
                boxfill(surf, base_x + 5, bar_y, p0_w, bar_h, rgb(surf->format, 80, 180, 100));
                boxfill(surf, base_x + 5 + p0_w, bar_y, bar_w - p0_w, bar_h, rgb(surf->format, 180, 80, 100));

                // Z expectation bar (simpler visualization)
                int z_bar_y = bar_y + 12;
                boxfill(surf, base_x + 5, z_bar_y, bar_w, 4, rgb(surf->format, 40, 40, 50));
                // Z ranges from -1 (|1>) to +1 (|0>)
                double z_norm = (v.Z + 1.0) / 2.0; // map to [0,1]
                int z_pos = base_x + 5 + (int)(z_norm * bar_w);
                boxfill(surf, z_pos - 2, z_bar_y - 2, 4, 8, rgb(surf->format, 230, 230, 240));
            }
        }

        SDL_UpdateWindowSurface(sdlm::M().bloch.win);

        Uint32 t0 = SDL_GetTicks();
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) { if (ev.type == SDL_QUIT) return false; }
        Uint32 dt = SDL_GetTicks() - t0;
        if ((int)dt < cfg.min_frame_ms) SDL_Delay(cfg.min_frame_ms - dt);
        return true;
    }

} // namespace quantum_viz