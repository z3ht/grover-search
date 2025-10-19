#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <cstddef>

namespace qcommon {

using Complex = std::complex<double>;
constexpr double PI = 3.1415926535897932384626433832795;

struct QRegister {
    int n = 0;                          // number of qubits
    std::vector<Complex> psi;           // statevector of length 2^n
    explicit QRegister(int n_bits)
            : n(n_bits), psi(1 << n_bits, Complex{0.0, 0.0}) {}
};

inline void hadamard_on_qubit(std::vector<Complex>& psi, int q) {
    const int N = static_cast<int>(psi.size());
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

inline void hadamard_all(std::vector<Complex>& psi, int n) {
    for (int q = 0; q < n; ++q) hadamard_on_qubit(psi, q);
}

inline void pauli_x_on_qubit(std::vector<Complex>& psi, int q) {
    const int N = static_cast<int>(psi.size());
    const int mask = 1 << q;
    for (int i = 0; i < N; ++i) {
        int j = i ^ mask;
        if (i < j) std::swap(psi[i], psi[j]);
    }
}

template <class Pred>
inline void oracle_phase_flip(std::vector<Complex>& psi, Pred pred) {
    const int N = static_cast<int>(psi.size());
    for (int i = 0; i < N; ++i) if (pred(i)) psi[i] = -psi[i];
}

inline void mcz_all_ones(std::vector<Complex>& psi, int n) {
    const int idx = (1 << n) - 1; // |11..1>
    psi[idx] = -psi[idx];
}

inline void diffusion_about_mean(std::vector<Complex>& psi) {
    qcommon::Complex mean{0.0, 0.0};
    for (const auto& a : psi) mean += a;
    mean /= static_cast<double>(psi.size());
    for (auto& a : psi) a = 2.0 * mean - a;
}

inline int grover_iterations_for(int N, int M) {
    if (M <= 0) return 0;
    double k = std::floor((PI / 4.0) * std::sqrt(static_cast<double>(N) / static_cast<double>(M)));
    return static_cast<int>(std::max(0.0, k));
}
inline double grover_theta_for(int N, int M) {
    return std::asin(std::sqrt(static_cast<double>(M) / static_cast<double>(N)));
}
inline double expected_marked_prob(int j, double theta) {
    const double ang = (2 * j + 1) * theta;
    return std::sin(ang) * std::sin(ang);
}

inline double l2_from_uniform(const std::vector<Complex>& psi) {
    const double N = static_cast<double>(psi.size());
    const qcommon::Complex u{1.0 / std::sqrt(N), 0.0};
    double acc = 0.0;
    for (const auto& a : psi) { const qcommon::Complex d = a - u; acc += std::norm(d); }
    return std::sqrt(acc);
}

inline bool bit_is_zero_in_target(int target, int q) { return ((target >> q) & 1) == 0; }

struct OracleMap {
    std::vector<int> flip_bits;   // which qubits to X
    int transformed_target = -1;  // target after pre-X
    int all_ones_index = -1;      // (1<<n)-1
};
inline OracleMap build_oracle_map(int n, int target) {
    OracleMap m;
    for (int q = 0; q < n; ++q) if (bit_is_zero_in_target(target, q)) m.flip_bits.push_back(q);
    m.transformed_target = target;
    for (int q : m.flip_bits) m.transformed_target ^= (1 << q);
    m.all_ones_index = (1 << n) - 1;
    return m;
}

}
