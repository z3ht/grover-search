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

namespace sdlm {
    struct Win { SDL_Window* win=nullptr; SDL_Surface* surf=nullptr; int W=0,H=0; };
    struct Manager { bool inited=false; Win grid, arrows, bars; };
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
            w.surf = SDL_GetWindowSurface(w.win); w.W=W; w.H=H; return true;
        }
        if (w.W != W || w.H != H) {
            SDL_DestroyWindow(w.win);
            w.win = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, W, H, SDL_WINDOW_SHOWN);
            if (!w.win) return false;
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
    ((Uint32*)((Uint8*)s->pixels + y*s->pitch))[x] = c;
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
inline void arrow(SDL_Surface* s, int x0,int y0,int x1,int y1, Uint32 shaft, Uint32 head,
                  double head_len=10.0, double head_ang_deg=23.0)
{
    line(s,x0,y0,x1,y1,shaft);
    double vx = x1-x0, vy = y1-y0, L = std::hypot(vx,vy);
    if (L < 1.0) return;
    double ux=vx/L, uy=vy/L, a = head_ang_deg * M_PI/180.0;
    double rx1 =  ux*std::cos(a) -  uy*std::sin(a);
    double ry1 =  ux*std::sin(a) +  uy*std::cos(a);
    double rx2 =  ux*std::cos(-a) - uy*std::sin(-a);
    double ry2 =  ux*std::sin(-a)+ uy*std::cos(-a);
    line(s, x1,y1, (int)std::round(x1 - rx1*head_len), (int)std::round(y1 - ry1*head_len), head);
    line(s, x1,y1, (int)std::round(x1 - rx2*head_len), (int)std::round(y1 - ry2*head_len), head);
}
inline void hsv2rgb(double h, double s, double v, Uint8& R, Uint8& G, Uint8& B){
    double r=0,g=0,b=0;
    if (s<=0){ r=g=b=v; }
    else {
        double hh = (h - std::floor(h)) * 6.0;
        int i = (int)std::floor(hh) % 6;
        double f = hh - std::floor(hh);
        double p = v*(1.0-s), q = v*(1.0-s*f), t = v*(1.0 - s*(1.0 - f));
        switch(i){
            case 0: r=v; g=t; b=p; break; case 1: r=q; g=v; b=p; break;
            case 2: r=p; g=v; b=t; break; case 3: r=p; g=q; b=v; break;
            case 4: r=t; g=p; b=v; break; case 5: r=v; g=p; b=q; break;
        }
    }
    R=(Uint8)std::round(std::clamp(r,0.0,1.0)*255.0);
    G=(Uint8)std::round(std::clamp(g,0.0,1.0)*255.0);
    B=(Uint8)std::round(std::clamp(b,0.0,1.0)*255.0);
}

struct GridArrowsConfig {
    int  r_bits = 0;
    int  cell   = 24;
    bool fill_prob = true;
    bool hue_by_phase = true;
    bool two_tone_sign = false;
    double gamma = 0.55;
    int  marked_index = -1;
    int  min_frame_ms = 33;
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
    SDL_FillRect(surf, nullptr, rgb(surf->format, 8,8,10));

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
            Uint8 R=0,G=0,B=0;
            if (cfg.hue_by_phase && !cfg.two_tone_sign){
                double hue = (th + M_PI) / (2.0*M_PI);
                hsv2rgb(hue, 1.0, v, R,G,B);
            } else {
                bool pos = (a.real() >= 0.0);
                if (pos){ R=(Uint8)(30*v*255); G=(Uint8)(90*v*255); B=(Uint8)(200*v*255); }
                else    { R=(Uint8)(220*v*255); G=(Uint8)(130*v*255); B=(Uint8)(35*v*255); }
            }
            SDL_Rect rc{ x0+1, y0+1, cfg.cell-2, cfg.cell-2 };
            SDL_FillRect(surf, &rc, rgb(surf->format, R,G,B));
        }

        double len01 = (cfg.arrow_len_mode == 0)
                       ? std::sqrt(std::clamp(amp,0.0,1.0))
                       : std::clamp(amp / max_amp, 0.0, 1.0);

        int tipx = (int)std::round(cx + std::cos(th) * len01 * (rad-2));
        int tipy = (int)std::round(cy - std::sin(th) * len01 * (rad-2));

        Uint8 rC=220,gC=220,bC=220;
        if (cfg.hue_by_phase && !cfg.two_tone_sign){
            double hue = (th + M_PI) / (2.0*M_PI);
            hsv2rgb(hue, 1.0, 1.0, rC,gC,bC);
        } else {
            if (a.real() >= 0.0){ rC=60; gC=160; bC=255; }
            else                { rC=255; gC=140; bC=50; }
        }
        Uint32 shaft = rgb(surf->format, rC,gC,bC);
        Uint32 head  = rgb(surf->format, std::min(255,(int)rC+12),
                           std::min(255,(int)gC+12),
                           std::min(255,(int)bC+12));
        arrow(surf, cx,cy, tipx,tipy, shaft, head);

        if (i == cfg.marked_index){
            Uint32 hi = rgb(surf->format, 255,60,60);
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
    bool   draw_theta_arc  = true;
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
inline void dashed_arc(SDL_Surface* s, int cx,int cy,int R, double ang0, double ang1, Uint32 col){
    const double step = 2.0 * M_PI / 360.0;
    double a = ang0;
    int on = 0;
    auto advance = [&](double& x){ x += (ang1 >= ang0 ? step : -step); };
    while ((ang1 >= ang0 && a <= ang1) || (ang1 < ang0 && a >= ang1)) {
        int x = (int)std::round(cx + std::cos(a) * R);
        int y = (int)std::round(cy - std::sin(a) * R);
        if ((on++ % 3) != 0) pset(s, x, y, col);
        advance(a);
    }
}

inline bool show_arrows(const std::vector<Complex>& psi, int n,
                        const ArrowsConfig& cfg, const char* title,
                        const std::vector<Complex>* prev_psi = nullptr)
{
    const int N = 1 << n;
    const int W = cfg.panel_size;
    const int H = cfg.panel_size + 48;

    if (!sdlm::ensure_window(sdlm::M().arrows, W, H, title)) return false;
    SDL_Surface* surf = sdlm::M().arrows.surf;
    SDL_FillRect(surf, nullptr, rgb(surf->format, 12,12,14));

    const int cx = W/2;
    const int cy = cfg.panel_size/2;
    const int R  = (int)std::round((cfg.panel_size/2 - cfg.margin) * cfg.arrow_scale);

    Uint32 axis = rgb(surf->format, 160,160,160);
    line(surf, cfg.margin, cy, cfg.panel_size - cfg.margin, cy, axis);
    line(surf, cx, cfg.margin, cx, cfg.panel_size - cfg.margin, axis);
    for (int deg=0; deg<360; ++deg){
        double t = deg * M_PI/180.0;
        pset(surf, (int)std::round(cx + R*std::cos(t)),
             (int)std::round(cy - R*std::sin(t)),
             rgb(surf->format, 70,70,90));
    }

    Complex mean(0,0);
    for (auto& a : psi) mean += a;
    mean /= (double)N;

    if (cfg.draw_all){
        Uint32 faint = rgb(surf->format, (int)(50*cfg.tail_alpha),
                           (int)(110*cfg.tail_alpha),
                           (int)(180*cfg.tail_alpha));
        for (int i=0;i<N;++i){
            if (i == cfg.marked_index) continue;
            auto [x1,y1] = to_xy(cx,cy,R,psi[i]);
            arrow(surf, cx,cy, x1,y1, faint, faint, 8.0, 20.0);
        }
    }

    if (cfg.draw_mean){
        auto [mx,my] = to_xy(cx,cy,R,mean);
        Uint32 g = rgb(surf->format, 40,200,120);
        arrow(surf, cx,cy, mx,my, g,g, 10.0, 22.0);
    }

    Complex pm = (cfg.marked_index>=0 && cfg.marked_index<N && prev_psi)
                 ? (*prev_psi)[cfg.marked_index] : Complex(0,0);
    if (prev_psi && cfg.marked_index >= 0 && cfg.marked_index < N){
        auto [px,py] = to_xy(cx,cy,R,pm);
        Uint32 red = rgb(surf->format, 230,80,70);
        arrow(surf, cx,cy, px,py, red, red, 8.0, 22.0);
    }

    if (cfg.marked_index >= 0 && cfg.marked_index < N){
        Complex am = psi[cfg.marked_index];
        auto [ax,ay] = to_xy(cx,cy,R,am);
        Uint32 cyan = rgb(surf->format, 60,200,255);
        arrow(surf, cx,cy, ax,ay, cyan, cyan, 12.0, 25.0);

        if (cfg.draw_reflection){
            Complex ref_src = prev_psi ? pm : am;
            Complex refl = 2.0*mean - ref_src;
            auto [rx,ry] = to_xy(cx,cy,R,refl);
            Uint32 magenta = rgb(surf->format, 220,100,240);
            int x0=cx, y0=cy, x1=rx, y1=ry;
            int dx=std::abs(x1-x0), sx=x0<x1?1:-1;
            int dy=-std::abs(y1-y0), sy=y0<y1?1:-1;
            int err = dx+dy, e2, step=0;
            while(true){
                if ((step++ % 3)!=0) pset(surf,x0,y0,magenta);
                if (x0==x1 && y0==y1) break;
                e2=2*err;
                if (e2>=dy){ err+=dy; x0+=sx; }
                if (e2<=dx){ err+=dx; y0+=sy; }
            }
        }

        if (cfg.draw_theta_arc){
            double ang_mean = std::arg(mean);
            double ang_mark = std::arg(psi[cfg.marked_index]);
            auto norm = [](double a){ while (a >= M_PI) a -= 2*M_PI; while (a < -M_PI) a += 2*M_PI; return a; };
            double d = norm(ang_mark - ang_mean);
            Uint32 col = rgb(surf->format, 200,200,60);
            dashed_arc(surf, cx, cy, (int)(R*0.75), ang_mean, ang_mean + d, col);
        }
    }

    int yb = cfg.panel_size + 8;
    auto box = [&](int x,int y, int r,int g,int b){
        SDL_Rect Rct{ x,y, 20,10 };
        SDL_FillRect(surf, &Rct, rgb(surf->format,r,g,b));
    };
    box(12,  yb,  60,200,255);
    box(112, yb,  40,200,120);
    box(212, yb, 230, 80, 70);
    box(332, yb, 220,100,240);
    box(452, yb, 200,200, 60);

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
    int    H = 320;
    bool   log_scale = false;
    double gamma = 0.6;
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
    SDL_FillRect(surf, nullptr, rgb(surf->format, 10,10,12));

    std::vector<std::pair<int,double>> items;
    items.reserve(N);
    for (int i=0;i<N;++i) items.emplace_back(i, std::norm(psi[i]));
    if (cfg.sort_descending) {
        std::stable_sort(items.begin(), items.end(),
                         [](auto& a, auto& b){ return a.second > b.second; });
    }

    Uint32 axis = rgb(surf->format, 140,140,140);
    line(surf, 32, cfg.H-28, cfg.W-16, cfg.H-28, axis);
    line(surf, 32, 16, 32, cfg.H-28, axis);

    const int left = 40, right = cfg.W - 16, bottom = cfg.H - 28, top = 16;
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

        Uint32 col;
        if (idx == cfg.marked_index) col = rgb(surf->format, 255, 90, 90);
        else {
            int r = (int)(60 + 120*v), g = (int)(120 + 100*v), b = (int)(255 * (0.35 + 0.65*v));
            col = rgb(surf->format, r,g,b);
        }

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

inline void shutdown(){
    sdlm::destroy_window(sdlm::M().grid);
    sdlm::destroy_window(sdlm::M().arrows);
    sdlm::destroy_window(sdlm::M().bars);
    if (sdlm::M().inited) SDL_Quit();
    sdlm::M().inited = false;
}

}
