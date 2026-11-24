#pragma once
#include "SurfGBS.hpp"

namespace libcdgbs {

    class Example {
    public:
        SurfGBS gbs;
        /// Reads the input, builds the surface and writes the result files.
        /// @p filename may carry a supported extension (.mgbs, .gbs, .cgb,
        /// .ngbs, .gbp, .mlp), which selects the reader; without one, the
        /// extensions are tried in that order. Returns false if no input
        /// could be read.
        bool say_hello(const std::string& filename = "", const SurfGBS::InputParams& params = SurfGBS::InputParams());
    };
}
