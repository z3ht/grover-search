#pragma once
#include <vector>
#include <string>
#include <random>
#include <numeric>     // std::gcd
#include <optional>
#include <tuple>
#include <complex>
#include <cmath>
#include <algorithm>
#include "ialgorithm.hpp"
#include "qcommon.hpp"
#include "frames.hpp"

struct ShorAlgorithm : public IAlgorithm {
    // Composite number to factor (small N like 15, 21, 33, 35, 91…)
    int N = 15;

    explicit ShorAlgorithm(int N_) : N(N_) {}

    // ---------- math helpers ----------
    static long long modmul(long long a, long long b, long long m) {
        return (a * b) % m; // fine for small demo N
    }
    static long long modpow(long long a, long long e, long long m) {
        long long r = 1 % m, base = a % m;
        while (e > 0) {
            if (e & 1) r = modmul(r, base, m);
            base = modmul(base, base, m);
            e >>= 1;
        }
        return r;
    }
    static int order_of_a_mod_N(int a, int N) {
        if (std::gcd(a, N) != 1) return -1;
        int r = 1;
        long long val = modpow(a, r, N);
        while (val != 1) {
            ++r;
            if (r > 4 * N) return -1; // guard for tiny demo
            val = modpow(a, r, N);
        }
        return r;
    }
    static std::pair<int, int> pick_a(int N) {
        std::mt19937 rng{std::random_device{}()};
        std::uniform_int_distribution<int> dist(2, N - 2);
        for (int tries = 0; tries < 64; ++tries) {
            int a = dist(rng);
            int g = std::gcd(a, N);
            if (g == 1) return {a, 1};
            if (g > 1 && g < N) return {a, g}; // immediate factor
        }
        for (int a = 2; a <= N - 2; ++a) {
            int g = std::gcd(a, N);
            if (g == 1) return {a, 1};
            if (g > 1 && g < N) return {a, g};
        }
        return {2, 1};
    }

    // ---------- viz helpers (pure state ops) ----------
    static void normalize(std::vector<qcommon::Complex>& psi) {
        double s = 0.0;
        for (auto &z : psi) s += std::norm(z);
        if (s == 0.0) return;
        const double scale = 1.0 / std::sqrt(s);
        for (auto &z : psi) z *= scale;
    }

    // Put register in uniform superposition |s>
    static void make_uniform(std::vector<qcommon::Complex>& psi, int n) {
        std::fill(psi.begin(), psi.end(), qcommon::Complex{0.0, 0.0});
        psi[0] = qcommon::Complex{1.0, 0.0};
        qcommon::hadamard_all(psi, n);
    }

    // Phase-encode f(x)=a^x mod N as angle 2π f(x)/N, with uniform magnitude
    static void encode_phase_modexp(std::vector<qcommon::Complex>& psi, int n, int a, int N) {
        const int D = 1 << n;
        const double inv = 1.0 / std::sqrt(static_cast<double>(D));
        for (int x = 0; x < D; ++x) {
            int fx = static_cast<int>(modpow(a, x, N));        // f(x) ∈ [0, N-1]
            double phi = 2.0 * qcommon::PI * (static_cast<double>(fx) / static_cast<double>(N));
            psi[x] = std::polar(inv, phi);
        }
    }

    // Mock readout distribution: place probability peaks at multiples of D/r
    // This mimics the QFT + measurement peaks you’d expect from periodicity.
    static void make_period_peaks(std::vector<qcommon::Complex>& psi, int n, int r) {
        const int D = 1 << n;
        std::fill(psi.begin(), psi.end(), qcommon::Complex{0.0, 0.0});
        if (r <= 0) { psi[0] = qcommon::Complex{1.0, 0.0}; return; }

        // Choose r bins at k*D/r (rounded) and assign equal amplitude there
        for (int m = 0; m < r; ++m) {
            int idx = static_cast<int>(std::llround( (static_cast<double>(m) * D) / r ));
            if (idx >= D) idx = D - 1;
            psi[idx] = qcommon::Complex{1.0, 0.0};
        }
        normalize(psi); // makes total probability 1
    }

    // ---------- IAlgorithm ----------
    std::vector<qframes::Frame> run_and_make_frames(qcommon::QRegister reg) override {
        std::vector<qframes::Frame> frames;

        const int n = reg.n;                  // keep caller’s n so viz matches
        const int D = 1 << n;
        reg.psi.assign(D, qcommon::Complex{0.0, 0.0});
        reg.psi[0] = qcommon::Complex{1.0, 0.0};
        int step = 0;

        // Pick a
        auto [a, gcd_hint] = pick_a(N);

        // Frame A — uniform superposition
        {
            auto psi = reg.psi;
            make_uniform(psi, n);
            frames.push_back(qframes::Frame{
                    psi, {},
                    "Shor: init — uniform superposition |s⟩",
                    "All basis states equally weighted",
                    "Equal probabilities",
                    "Algorithm",
                    false, false, qframes::PhaseTag::None, step++, -1
            });
        }

        // Early luck: gcd gives factor immediately
        if (gcd_hint != 1 && gcd_hint != 0 && gcd_hint != N) {
            int p = gcd_hint, q = N / gcd_hint;
            frames.push_back(qframes::Frame{
                    reg.psi, {},
                    "Shor: lucky — gcd(a,N) > 1",
                    "a = " + std::to_string(a) + ", gcd(a,N) = " + std::to_string(gcd_hint),
                    "Factors: " + std::to_string(p) + " × " + std::to_string(q),
                    "Algorithm",
                    false, false, qframes::PhaseTag::None, step++, -1
            });
            return frames;
        }

        // Frame B — phase-encode f(x)=a^x mod N
        {
            auto psi = std::vector<qcommon::Complex>(D);
            encode_phase_modexp(psi, n, a, N);
            frames.push_back(qframes::Frame{
                    psi, {},
                    "Shor: phase-encode  f(x)=a^x mod N",
                    "a = " + std::to_string(a) + ", N = " + std::to_string(N),
                    "Phases carry periodicity (Argand shows spiral)",
                    "Algorithm",
                    false, false, qframes::PhaseTag::None, step++, -1
            });
        }

        // Find period r (demo: classically) and show a mock readout with r peaks
        int r = order_of_a_mod_N(a, N);
        {
            auto psi = std::vector<qcommon::Complex>(D);
            make_period_peaks(psi, n, r);
            frames.push_back(qframes::Frame{
                    psi, {},
                    "Shor: (mock) readout pattern",
                    "Estimated period r = " + std::to_string(r),
                    "Peaks at multiples of (2^n)/r",
                    "Algorithm",
                    false, false, qframes::PhaseTag::None, step++, -1
            });
        }

        auto factors = try_factor_from_period(a, N, r);
        if (factors) {
            frames.push_back(qframes::Frame{
                    reg.psi, {},
                    "Shor: factors",
                    "p = " + std::to_string(factors->first) + ", q = " + std::to_string(factors->second),
                    "Done",
                    "Algorithm",
                    false, false, qframes::PhaseTag::None, step++, -1
            });
        } else {
            frames.push_back(qframes::Frame{
                    reg.psi, {},
                    "Shor: note",
                    "If r is odd or a^(r/2) ≡ ±1 (mod N), pick a new a and retry.",
                    "Try different a for real run",
                    "Algorithm",
                    false, false, qframes::PhaseTag::None, step++, -1
            });
        }

        return frames;
    }

private:
    static std::optional<std::pair<int,int>> try_factor_from_period(int a, int N, int r) {
        if (r <= 0 || (r % 2) == 1) return std::nullopt;
        long long ar2 = modpow(a, r / 2, N);
        long long modv = ar2 % N;
        if (modv == 1 || modv == N - 1) return std::nullopt;
        int p = std::gcd(static_cast<int>(ar2 + 1), N);
        int q = std::gcd(static_cast<int>(ar2 - 1), N);
        if (p == 1 || p == N || q == 1 || q == N) return std::nullopt;
        return std::make_pair(p, q);
    }
};
