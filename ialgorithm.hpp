#pragma once

#include <vector>
#include "qcommon.hpp"
#include "frames.hpp"


struct IAlgorithm {
    virtual ~IAlgorithm() = default;
    virtual std::vector<qframes::Frame> run_and_make_frames(qcommon::QRegister initial) = 0;
};
