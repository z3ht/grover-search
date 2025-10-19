#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <cmath>
#include <bitset>
#include <algorithm>


using Complex = std::complex<double>;
constexpr double PI = 3.1415926535897932384626433832795;

void print_state(const std::vector<Complex>& psi, const std::string& label, double threshold = 1e-12) {
    std::cout << "\n=== " << label << " ===\n";
    int n = 0;
    while ((1 << n) < static_cast<int>(psi.size())) ++n;

    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(6);

    for (int i = 0; i < static_cast<int>(psi.size()); ++i) {
        double p = std::norm(psi[i]);
        if (p > threshold) {
            std::cout << std::setw(2) << i << " |"
                      << std::bitset<20>(i).to_string().substr(20 - n)
                      << ">: amp=" << psi[i].real()
                      << (psi[i].imag() >= 0 ? "+" : "")
                      << psi[i].imag() << "i, P=" << p << "\n";
        }
    }
}

double state_diff(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        Complex d = a[i] - b[i];
        s += std::norm(d);
    }
    return std::sqrt(s);
}

// ---------------------------------------------
// Quantum gates (state-vector manipulation)
// ---------------------------------------------
void apply_hadamard_on_qubit(std::vector<Complex>& psi, int q, int n) {
    const int N = static_cast<int>(psi.size());
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    const int stride = 1 << (q + 1);
    const int half = 1 << q;

    for (int base = 0; base < N; base += stride) {
        for (int i = 0; i < half; ++i) {
            int i0 = base + i;
            int i1 = base + i + half;
            Complex a = psi[i0];
            Complex b = psi[i1];
            psi[i0] = (a + b) * inv_sqrt2;
            psi[i1] = (a - b) * inv_sqrt2;
        }
    }
}

void apply_hadamard_all(std::vector<Complex>& psi, int n) {
    for (int q = 0; q < n; ++q)
        apply_hadamard_on_qubit(psi, q, n);
}

void apply_pauli_x_on_qubit(std::vector<Complex>& psi, int q, int n) {
    const int N = static_cast<int>(psi.size());
    const int mask = 1 << q;
    for (int i = 0; i < N; ++i) {
        int j = i ^ mask;
        if (i < j)
            std::swap(psi[i], psi[j]);
    }
}

void apply_x_all(std::vector<Complex>& psi, int n) {
    for (int q = 0; q < n; ++q)
        apply_pauli_x_on_qubit(psi, q, n);
}

// Multi-controlled Z on |11..1>
void apply_mcz_all_ones(std::vector<Complex>& psi, int n) {
    int idx = (1 << n) - 1;
    psi[idx] *= -1.0;
}

// ---------------------------------------------
// Grover primitives
// ---------------------------------------------
template <class Pred>
void apply_oracle_phase_flip(std::vector<Complex>& psi, Pred f) {
    const int N = static_cast<int>(psi.size());
    for (int i = 0; i < N; ++i)
        if (f(i)) psi[i] = -psi[i];
}

// Diffusion: algebraic inversion about the mean
void apply_diffusion_mean(std::vector<Complex>& psi) {
    Complex mean = {0.0, 0.0};
    for (const auto& a : psi) mean += a;
    mean /= static_cast<double>(psi.size());
    for (auto& a : psi) a = 2.0 * mean - a;
}

// Diffusion via gate decomposition
void apply_diffusion_gates(std::vector<Complex>& psi, int n) {
    apply_hadamard_all(psi, n);
    apply_x_all(psi, n);
    apply_mcz_all_ones(psi, n);
    apply_x_all(psi, n);
    apply_hadamard_all(psi, n);
}

// Recommended number of iterations
int grover_iterations(int N, int M) {
    if (M <= 0) return 0;
    double k = std::floor((PI / 4.0) * std::sqrt(static_cast<double>(N) / static_cast<double>(M)));
    return static_cast<int>(std::max(0.0, k));
}

// ---------------------------------------------
// Main simulation
// ---------------------------------------------
int main() {
    const int n = 3;
    const int N = 1 << n;

    auto is_marked = [](int x) -> bool {
        return (x == 5);
    };

    int M = 0;
    for (int i = 0; i < N; ++i)
        if (is_marked(i)) ++M;

    if (M == 0) {
        std::cerr << "No marked items; exiting.\n";
        return 0;
    }

    std::vector<Complex> psi(N, Complex{0.0, 0.0});
    psi[0] = Complex{1.0, 0.0};
    apply_hadamard_all(psi, n);
    print_state(psi, "After H^n (uniform |s>)");

    const int k = grover_iterations(N, M);
    std::cout << "\nN=" << N << ", M=" << M
              << ", recommended iterations k=" << k << "\n";

    const bool use_mean_diffusion = true;
    std::vector<Complex> psi_ref;

    for (int iter = 1; iter <= k; ++iter) {
        apply_oracle_phase_flip(psi, is_marked);
        print_state(psi, "After Oracle U_omega, iter " + std::to_string(iter));

        if (use_mean_diffusion) {
            psi_ref = psi;
            apply_diffusion_mean(psi);

            std::vector<Complex> psi_gate = psi_ref;
            apply_diffusion_gates(psi_gate, n);
            double diff = state_diff(psi, psi_gate);
            std::cout << "Diffusion mean vs gates L2 diff: "
                      << std::setprecision(3) << std::scientific << diff << "\n";
        } else {
            apply_diffusion_gates(psi, n);
        }

        print_state(psi, "After Diffusion D, iter " + std::to_string(iter));
    }

    print_state(psi, "Final state before measurement");

    std::cout << "\nTop probabilities:\n";
    std::vector<std::pair<double, int>> probs;
    probs.reserve(N);
    for (int i = 0; i < N; ++i)
        probs.emplace_back(std::norm(psi[i]), i);

    std::sort(probs.begin(), probs.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    for (int i = 0; i < std::min(N, 8); ++i) {
        std::cout << "x=" << std::setw(2) << probs[i].second
                  << "  P=" << std::fixed << std::setprecision(6)
                  << probs[i].first
                  << (is_marked(probs[i].second) ? "   <-- MARKED\n" : "\n");
    }
}
