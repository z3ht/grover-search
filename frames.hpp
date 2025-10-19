#pragma once
#include <vector>
#include <string>
#include "qcommon.hpp"


namespace qframes {

enum class PhaseTag : int {
    None = 0,
    PreX,
    MCZ,
    UnX,
    OracleResult,
};

struct Frame {
    std::vector<qcommon::Complex> psi;
    std::vector<qcommon::Complex> prev_for_reflection;
    std::string title_grid, title_arrows, title_bars, title_circuit;
    bool draw_reflection = false;
    bool show_circuit    = false;
    PhaseTag phase       = PhaseTag::None;
    int  iter_index      = 0;
    int  highlight_index = -1;
};

}
