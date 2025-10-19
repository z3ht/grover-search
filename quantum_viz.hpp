#pragma once
#include <SDL.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace quantum_viz {

    using Complex = std::complex<double>;
    constexpr double PI = 3.1415926535897932384626433832795;

    namespace sdlm {
        struct Win { SDL_Window* win=nullptr; SDL_Surface* surf=nullptr; int W=0,H=0; };
        struct Manager { bool inited=false; Win grid, arrows, bars, circuit, bloch; };
        inline Manager& M(){ static Manager m; return m; }

        inline bool ensure_sdl() {
            auto& m = M();
            if (!m.inited) { if (SDL_Init(SDL_INIT_VIDEO) < 0) return false; m.inited = true; }
            return true;
        }
        inline bool ensure_window(Win& w, int W, int H, const char* title) {
            if (!ensure_sdl()) return false;
            if (!w.win) {
                w.win = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, W, H, SDL_WINDOW_SHOWN);
                if (!w.win) return false;
                SDL_SetWindowResizable(w.win, SDL_FALSE); // non-resizable to avoid corner artifact
                w.surf = SDL_GetWindowSurface(w.win); w.W=W; w.H=H; return true;
            }
            if (w.W != W || w.H != H) {
                SDL_DestroyWindow(w.win);
                w.win = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, W, H, SDL_WINDOW_SHOWN);
                if (!w.win) return false;
                SDL_SetWindowResizable(w.win, SDL_FALSE);
                w.surf = SDL_GetWindowSurface(w.win); w.W=W; w.H=H;
            }
            SDL_SetWindowTitle(w.win, title);
            return true;
        }
        inline void destroy_window(Win& w) {
            if (w.win) SDL_DestroyWindow(w.win);
            w.win = nullptr; w.surf = nullptr; w.W = w.H = 0;
        }
    }

    inline Uint32 rgb(SDL_PixelFormat* f, int r,int g,int b) {
        return SDL_MapRGB(f, (Uint8)std::clamp(r,0,255), (Uint8)std::clamp(g,0,255), (Uint8)std::clamp(b,0,255));
    }
    inline void pset(SDL_Surface* s, int x,int y, Uint32 c){
        if (!s) return;
        if (x<0||y<0||x>=s->w||y>=s->h) return;
        reinterpret_cast<Uint32*>(static_cast<Uint8*>(s->pixels) + y * s->pitch)[x] = c;
    }
    inline void line(SDL_Surface* s, int x0,int y0,int x1,int y1, Uint32 c){
        int dx=std::abs(x1-x0), sx=x0<x1?1:-1;
        int dy=-std::abs(y1-y0), sy=y0<y1?1:-1;
        int err=dx+dy, e2;
        for(;;){
            pset(s,x0,y0,c);
            if (x0==x1 && y0==y1) break;
            e2=2*err;
            if (e2>=dy){ err+=dy; x0+=sx; }
            if (e2<=dx){ err+=dx; y0+=sy; }
        }
    }
    inline void boxfill(SDL_Surface* s, int x,int y,int w,int h, Uint32 c){
        SDL_Rect rc{ x,y,w,h }; SDL_FillRect(s, &rc, c);
    }
    inline void circle(SDL_Surface* s, int cx,int cy,int R, Uint32 c){
        for (int a=0;a<360;++a){
            double t=a*PI/180.0;
            pset(s, (int)std::round(cx+R*std::cos(t)), (int)std::round(cy-R*std::sin(t)), c);
        }
    }
    inline void arrow(SDL_Surface* s, int x0,int y0,int x1,int y1, Uint32 shaft, Uint32 head,
                      double head_len=10.0, double head_ang_deg=23.0)
    {
        line(s,x0,y0,x1,y1,shaft);
        double vx = x1-x0, vy = y1-y0, L = std::hypot(vx,vy);
        if (L < 1.0) return;
        double ux=vx/L, uy=vy/L, a = head_ang_deg * PI/180.0;
        double rx1 =  ux*std::cos(a) -  uy*std::sin(a);
        double ry1 =  ux*std::sin(a) +  uy*std::cos(a);
        double rx2 =  ux*std::cos(-a) - uy*std::sin(-a);
        double ry2 =  ux*std::sin(-a)+ uy*std::cos(a);
        line(s, x1,y1, (int)std::round(x1 - rx1*head_len), (int)std::round(y1 - ry1*head_len), head);
        line(s, x1,y1, (int)std::round(x1 - rx2*head_len), (int)std::round(y1 - ry2*head_len), head);
    }

    struct GridArrowsConfig {
        int  r_bits = 0;
        int  cell   = 24;
        bool fill_prob = true;
        bool hue_by_phase = false;
        bool two_tone_sign = false;
        double gamma = 0.55;
        int  marked_index = -1;
        int  min_frame_ms = 16;
        double arrow_len_mode = 0;
    };

    inline bool show_grid_arrows(const std::vector<Complex>& psi, int n,
                                 const GridArrowsConfig& cfg, const char* title)
    {
        const int N = 1 << n;
        const int rows = 1 << cfg.r_bits;
        const int cols = 1 << (n - cfg.r_bits);
        const int W = cols * cfg.cell;
        const int H = rows * cfg.cell;

        if (!sdlm::ensure_window(sdlm::M().grid, W, H, title)) return false;
        SDL_Surface* surf = sdlm::M().grid.surf;
        SDL_FillRect(surf, nullptr, rgb(surf->format, 12,12,14));

        auto idx_to_rc = [&](int i, int& r, int& c){
            int c_bits = n - cfg.r_bits;
            r = (c_bits>0) ? (i >> c_bits) : i;
            int mask = (c_bits>0) ? ((1<<c_bits)-1) : 0;
            c = (c_bits>0) ? (i & mask) : 0;
        };

        double max_amp = 0.0;
        for (auto& a : psi) max_amp = std::max(max_amp, std::abs(a));
        if (max_amp <= 0.0) max_amp = 1.0;

        for (int i=0; i<N; ++i){
            int r,c; idx_to_rc(i,r,c);
            int x0 = c*cfg.cell, y0 = r*cfg.cell;
            int cx = x0 + cfg.cell/2;
            int cy = y0 + cfg.cell/2;
            int rad = (int)std::round(cfg.cell*0.5 - 2.0);

            const Complex a = psi[i];
            const double p  = std::norm(a);
            const double amp = std::abs(a);
            const double th  = std::atan2(a.imag(), a.real());

            if (cfg.fill_prob){
                double v = std::pow(std::clamp(p,0.0,1.0), cfg.gamma);
                Uint8 g = (Uint8)std::round(255.0 * v * 0.6);
                SDL_Rect rc{ x0+1, y0+1, cfg.cell-2, cfg.cell-2 };
                SDL_FillRect(surf, &rc, rgb(surf->format, g,g,g));
            }

            double len01 = (cfg.arrow_len_mode == 0)
                           ? std::sqrt(std::clamp(amp,0.0,1.0))
                           : std::clamp(amp / max_amp, 0.0, 1.0);

            int tipx = (int)std::round(cx + std::cos(th) * len01 * (rad-2));
            int tipy = (int)std::round(cy - std::sin(th) * len01 * (rad-2));

            bool is_marked = (i == cfg.marked_index);
            Uint32 shaft = is_marked ? rgb(surf->format, 50,200,255) : rgb(surf->format, 200,200,200);
            Uint32 head  = is_marked ? rgb(surf->format, 80,220,255) : rgb(surf->format, 220,220,220);
            arrow(surf, cx,cy, tipx,tipy, shaft, head);

            if (is_marked){
                Uint32 hi = rgb(surf->format, 50,200,255);
                for (int dx=0; dx<cfg.cell; ++dx){ pset(surf, x0+dx, y0, hi); pset(surf, x0+dx, y0+cfg.cell-1, hi); }
                for (int dy=0; dy<cfg.cell; ++dy){ pset(surf, x0, y0+dy, hi); pset(surf, x0+cfg.cell-1, y0+dy, hi); }
            }
        }

        SDL_UpdateWindowSurface(sdlm::M().grid.win);

        Uint32 t0 = SDL_GetTicks();
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) { if (ev.type == SDL_QUIT) return false; }
        Uint32 dt = SDL_GetTicks() - t0;
        if ((int)dt < cfg.min_frame_ms) SDL_Delay(cfg.min_frame_ms - dt);
        return true;
    }

    struct ArrowsConfig {
        int    panel_size   = 640;
        int    margin       = 16;
        int    min_frame_ms = 16;
        int    marked_index = -1;
        bool   draw_all     = true;
        bool   draw_mean    = true;
        bool   draw_reflection = false;
        bool   draw_theta_arc  = false;
        double tail_alpha   = 0.28;
        double arrow_scale  = 0.9;
    };

    inline std::pair<int,int> to_xy(int cx,int cy,int R, Complex a){
        double r = std::abs(a);
        double th = std::arg(a);
        int x = (int)std::round(cx + std::cos(th) * r * R);
        int y = (int)std::round(cy - std::sin(th) * r * R);
        return {x,y};
    }

    inline bool show_arrows(const std::vector<Complex>& psi, int n,
                            const ArrowsConfig& cfg, const char* title,
                            const std::vector<Complex>* prev_psi = nullptr)
    {
        const int N = 1 << n;
        const int W = cfg.panel_size;
        const int H = cfg.panel_size + 30;

        if (!sdlm::ensure_window(sdlm::M().arrows, W, H, title)) return false;
        SDL_Surface* surf = sdlm::M().arrows.surf;
        SDL_FillRect(surf, nullptr, rgb(surf->format, 12,12,14));

        const int cx = W/2;
        const int cy = cfg.panel_size/2;
        const int R  = (int)std::round((cfg.panel_size/2 - cfg.margin) * cfg.arrow_scale);

        Uint32 axis = rgb(surf->format, 180,180,180);
        line(surf, cfg.margin, cy, cfg.panel_size - cfg.margin, cy, axis);
        line(surf, cx, cfg.margin, cx, cfg.panel_size - cfg.margin, axis);
        circle(surf, cx, cy, R, rgb(surf->format, 60,60,70));

        Complex mean(0,0);
        for (auto& a : psi) mean += a;
        mean /= (double)N;

        if (cfg.draw_all){
            Uint32 faint = rgb(surf->format, (int)(200*cfg.tail_alpha),
                               (int)(200*cfg.tail_alpha),
                               (int)(200*cfg.tail_alpha));
            for (int i=0;i<N;++i){
                if (i == cfg.marked_index) continue;
                auto [x1,y1] = to_xy(cx,cy,R,psi[i]);
                arrow(surf, cx,cy, x1,y1, faint, faint, 8.0, 20.0);
            }
        }

        if (cfg.draw_mean){
            auto [mx,my] = to_xy(cx,cy,R,mean);
            Uint32 g = rgb(surf->format, 200,200,200);
            arrow(surf, cx,cy, mx,my, g,g, 10.0, 22.0);
        }

        Complex pm = (cfg.marked_index>=0 && cfg.marked_index<N && prev_psi)
                     ? (*prev_psi)[cfg.marked_index] : Complex(0,0);
        if (prev_psi && cfg.marked_index >= 0 && cfg.marked_index < N){
            auto [px,py] = to_xy(cx,cy,R,pm);
            Uint32 ghost = rgb(surf->format, 140,140,160);
            arrow(surf, cx,cy, px,py, ghost, ghost, 8.0, 22.0);
        }

        if (cfg.marked_index >= 0 && cfg.marked_index < N){
            Complex am = psi[cfg.marked_index];
            auto [ax,ay] = to_xy(cx,cy,R,am);
            Uint32 cyan = rgb(surf->format, 50,200,255);
            arrow(surf, cx,cy, ax,ay, cyan, cyan, 12.0, 25.0);

            if (cfg.draw_reflection){
                Complex ref_src = prev_psi ? pm : am;
                Complex refl = 2.0*mean - ref_src;
                auto [rx,ry] = to_xy(cx,cy,R,refl);
                Uint32 cyan_d = rgb(surf->format, 80,220,255);
                int x0=cx, y0=cy, x1=rx, y1=ry;
                int dx=std::abs(x1-x0), sx=x0<x1?1:-1;
                int dy=-std::abs(y1-y0), sy=y0<y1?1:-1;
                int err = dx+dy, e2, step=0;
                while(true){
                    if ((step++ % 3)!=0) pset(surf,x0,y0,cyan_d);
                    if (x0==x1 && y0==y1) break;
                    e2=2*err;
                    if (e2>=dy){ err+=dy; x0+=sx; }
                    if (e2<=dx){ err+=dx; y0+=sy; }
                }
            }
        }

        SDL_UpdateWindowSurface(sdlm::M().arrows.win);

        Uint32 t0 = SDL_GetTicks();
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) { if (ev.type == SDL_QUIT) return false; }
        Uint32 dt = SDL_GetTicks() - t0;
        if ((int)dt < cfg.min_frame_ms) SDL_Delay(cfg.min_frame_ms - dt);
        return true;
    }

    struct ProbBarsConfig {
        int    W = 640;
        int    H = 300;
        bool   log_scale = false;
        double gamma = 0.55;
        bool   sort_descending = false;
        int    marked_index = -1;
        int    min_frame_ms = 16;
    };

    inline bool show_prob_bars(const std::vector<Complex>& psi, int n,
                               const ProbBarsConfig& cfg, const char* title)
    {
        const int N = 1 << n;
        if (!sdlm::ensure_window(sdlm::M().bars, cfg.W, cfg.H, title)) return false;
        SDL_Surface* surf = sdlm::M().bars.surf;
        SDL_FillRect(surf, nullptr, rgb(surf->format, 12,12,14));

        std::vector<std::pair<int,double>> items;
        items.reserve(N);
        for (int i=0;i<N;++i) items.emplace_back(i, std::norm(psi[i]));
        if (cfg.sort_descending) {
            std::stable_sort(items.begin(), items.end(),
                             [](auto& a, auto& b){ return a.second > b.second; });
        }

        Uint32 axis = rgb(surf->format, 180,180,180);
        line(surf, 32, cfg.H-24, cfg.W-16, cfg.H-24, axis);
        line(surf, 32, 16, 32, cfg.H-24, axis);

        const int left = 40, right = cfg.W - 16, bottom = cfg.H - 24, top = 16;
        const int span = right - left;
        const double bar_w_f = std::max(1.0, span / (double)N);
        const int bar_w = (int)std::floor(bar_w_f);

        double maxp = 0.0;
        for (auto& kv: items) maxp = std::max(maxp, kv.second);
        if (maxp <= 0.0) maxp = 1.0;

        for (int k=0; k<N; ++k) {
            int idx = items[k].first;
            double p  = items[k].second;

            double v = cfg.log_scale ? (std::log1p(p) / std::log1p(maxp)) : (p / maxp);
            v = std::pow(std::clamp(v, 0.0, 1.0), cfg.gamma);

            int x0 = left + k*bar_w;
            int x1 = x0 + bar_w - 1;
            int h  = (int)std::round(v * (bottom - top));
            int y0 = bottom - h;

            Uint32 col = (idx == cfg.marked_index)
                         ? rgb(surf->format, 50,200,255)
                         : rgb(surf->format, 200,200,200);

            SDL_Rect rc{ x0, y0, std::max(1, x1-x0+1), h };
            SDL_FillRect(surf, &rc, col);
        }

        SDL_UpdateWindowSurface(sdlm::M().bars.win);

        Uint32 t0 = SDL_GetTicks();
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) { if (ev.type == SDL_QUIT) return false; }
        Uint32 dt = SDL_GetTicks() - t0;
        if ((int)dt < cfg.min_frame_ms) SDL_Delay(cfg.min_frame_ms - dt);
        return true;
    }

    struct BlochSpheresConfig {
        int sphere_size = 180;
        int margin = 20;
        int cols = 4;
        int min_frame_ms = 16;
        bool draw_axes = true;
        bool draw_state_vector = true;
        bool draw_probabilities = true;
    };

    bool show_bloch_spheres(const std::vector<Complex>& psi, int n,
                            const BlochSpheresConfig& cfg, const char* title);

    enum class Step { None, Prev, Next, First, Last, Quit };
    inline Step wait_for_step(){
        SDL_Event ev;
        for (;;) {
            while (SDL_PollEvent(&ev)) {
                if (ev.type == SDL_QUIT) return Step::Quit;
                if (ev.type == SDL_KEYDOWN) {
                    switch (ev.key.keysym.sym) {
                        case SDLK_LEFT:  case SDLK_a:    return Step::Prev;
                        case SDLK_RIGHT: case SDLK_d:    return Step::Next;
                        case SDLK_HOME:                  return Step::First;
                        case SDLK_END:                   return Step::Last;
                        case SDLK_ESCAPE: case SDLK_q:  return Step::Quit;
                        default: break;
                    }
                }
            }
            SDL_Delay(12);
        }
    }

    enum class OraclePhase { None, PreX, MCX, UnX, Result };

    struct OraclePanelConfig {
        int W = 720;
        int H = 240;
        int margin = 28;
        int min_frame_ms = 16;
    };

    inline bool show_oracle_step(int n, int target, OraclePhase phase,
                                 const OraclePanelConfig& cfg, const char* title)
    {
        if (!sdlm::ensure_window(sdlm::M().circuit, cfg.W, cfg.H, title)) return false;
        SDL_Surface* surf = sdlm::M().circuit.surf;
        SDL_FillRect(surf, nullptr, rgb(surf->format, 12,12,14));

        const int lines = n + 1;
        const int top = cfg.margin;
        const int left = cfg.margin;
        const int right = cfg.W - cfg.margin;
        const int vstep = (cfg.H - 2*cfg.margin) / std::max(1, lines - 1);

        Uint32 wire = rgb(surf->format, 180,180,180);
        Uint32 cyan = rgb(surf->format, 50,200,255);
        Uint32 gate = rgb(surf->format, 200,200,200);

        for (int i=0;i<lines;++i){
            int y = top + i*vstep;
            line(surf, left, y, right, y, wire);
        }

        auto draw_box = [&](int x,int y, Uint32 col){ boxfill(surf, x-14, y-9, 28, 18, col); };
        auto draw_dot = [&](int x,int y, Uint32 col){ boxfill(surf, x-3, y-3, 6, 6, col); };
        auto draw_plus = [&](int x,int y, Uint32 col){
            circle(surf, x, y, 8, col); line(surf, x-6,y, x+6,y, col); line(surf, x,y-6, x,y+6, col);
        };
        auto vbar = [&](int x,int y0,int y1, Uint32 col){ line(surf, x, y0, x, y1, col); };

        int col_x1 = left + (right-left)*2/6;
        int col_x2 = left + (right-left)*3/6;
        int col_x3 = left + (right-left)*4/6;

        int y_head = top - 8;
        auto head_tick = [&](int x, Uint32 col){ boxfill(surf, x-2, y_head-10, 4, 10, col); };
        head_tick(col_x1, cyan);
        head_tick(col_x2, (phase==OraclePhase::MCX)? cyan : wire);
        head_tick(col_x3, (phase==OraclePhase::UnX)? cyan : wire);

        int y_bits = top - 22;
        for (int q=0; q<n; ++q){
            bool bit1 = ((target >> q) & 1);
            int x = left + 10 + q * 12;
            Uint32 fill = bit1 ? rgb(surf->format,200,200,200) : rgb(surf->format,60,60,60);
            boxfill(surf, x-5, y_bits-5, 10, 10, fill);
            Uint32 ol = rgb(surf->format,120,120,120);
            line(surf, x-5, y_bits-5, x+5, y_bits-5, ol);
            line(surf, x-5, y_bits+5, x+5, y_bits+5, ol);
            line(surf, x-5, y_bits-5, x-5, y_bits+5, ol);
            line(surf, x+5, y_bits-5, x+5, y_bits+5, ol);
        }

        if (phase == OraclePhase::PreX || phase == OraclePhase::UnX){
            for (int q=0;q<n;++q){
                bool flip = (((target >> q) & 1) == 0);
                if (flip){
                    int y = top + q*vstep;
                    draw_box((phase==OraclePhase::PreX)? col_x1 : col_x3, y, cyan);
                }
            }
        }

        if (phase == OraclePhase::MCX){
            int y_lo = top;
            int y_hi = top + (n-1)*vstep;
            vbar(col_x2, y_lo, y_hi, wire);
            for (int q=0;q<n;++q){
                int y = top + q*vstep;
                draw_dot(col_x2, y, gate);
            }
            int y_anc = top + n*vstep;
            draw_plus(col_x2, y_anc, cyan);
        }

        SDL_UpdateWindowSurface(sdlm::M().circuit.win);

        Uint32 t0 = SDL_GetTicks();
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) { if (ev.type == SDL_QUIT) return false; }
        Uint32 dt = SDL_GetTicks() - t0;
        if ((int)dt < cfg.min_frame_ms) SDL_Delay(cfg.min_frame_ms - dt);
        return true;
    }

    inline void shutdown(){
        sdlm::destroy_window(sdlm::M().grid);
        sdlm::destroy_window(sdlm::M().arrows);
        sdlm::destroy_window(sdlm::M().bars);
        sdlm::destroy_window(sdlm::M().circuit);
        sdlm::destroy_window(sdlm::M().bloch);
        if (sdlm::M().inited) SDL_Quit();
        sdlm::M().inited = false;
    }

}
