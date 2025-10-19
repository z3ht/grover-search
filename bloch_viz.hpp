#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <cstdio>
#include <string>

namespace qviz_local
{
    struct QubitViz
    {
        double p0, p1;
        double X, Y, Z;
        double theta, phi;
        double length;
    };

    inline std::vector<QubitViz>
    compute_qubit_bloch_vectors(const std::vector<std::complex<double>>& psi, int n_qubits)
    {
        using C = std::complex<double>;
        const size_t N = psi.size();
        std::vector<QubitViz> out(n_qubits);

        for (int q = 0; q < n_qubits; ++q)
        {
            const size_t bit = (1ULL << q);
            double rho00 = 0.0, rho11 = 0.0;
            C rho01 = 0.0;

            for (size_t base = 0; base < N; ++base)
            {
                if (base & bit) continue;       // only where bit q = 0
                size_t i0 = base;
                size_t i1 = base | bit;
                const C a0 = psi[i0];
                const C a1 = psi[i1];

                rho00 += std::norm(a0);
                rho11 += std::norm(a1);
                rho01 += a0 * std::conj(a1);
            }

            const double X =  2.0 * std::real(rho01);
            const double Y = -2.0 * std::imag(rho01);
            const double Z =  rho00 - rho11;
            const double len = std::sqrt(std::max(0.0, X*X + Y*Y + Z*Z));
            const double theta = (len > 1e-12) ? std::acos(std::clamp(Z / len, -1.0, 1.0)) : 0.0;
            const double phi   = std::atan2(Y, X);

            out[q] = QubitViz{rho00, rho11, X, Y, Z, theta, phi, len};
        }

        return out;
    }

    inline std::string barZ(double Z, int width = 21)
    {
        int mid = width / 2;
        int pos = mid + int(std::round(Z * mid));
        if (pos < 0) pos = 0;
        if (pos >= width) pos = width - 1;
        std::string s(width, '.');
        s[pos] = '|';
        return s;
    }

    inline void print_qubit_summary(const std::vector<std::complex<double>>& psi, int n_qubits)
    {
        auto qubits = compute_qubit_bloch_vectors(psi, n_qubits);
        for (int q = 0; q < n_qubits; ++q)
        {
            const auto& v = qubits[q];
            std::printf("q%-2d  Z:%s  p0=%.3f p1=%.3f |B|=%.2f θ=%.2f φ=%.2f\n",
                        q, barZ(v.Z).c_str(), v.p0, v.p1, v.length, v.theta, v.phi);
        }
    }
}
